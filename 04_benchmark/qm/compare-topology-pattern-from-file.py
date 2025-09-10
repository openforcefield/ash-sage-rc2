"""
This script compares topology values between QM and MM for a single topology pattern.
This was constructed to run efficiently as a batch script on a SLURM cluster,
so it reads a file of topology patterns, and only runs for the pattern at the specified index..
\b
The input data directory is expected to contain parquet files
with topology values computed by the `compute-topology-values.py` script.
The input JSON file is expected to contain a list of dictionaries.
Each dictionary should contain the following keys:
    - pattern (str): The SMARTS pattern to search for.
    - output_stem (str): The stem for the output files.
    - topology_group (str): The topology group to search for, e.g., "Bonds", "Angles", "ImproperTorsions".

The output files are written using the specified `output_stem` from the input JSON file.
Two files are written: a CSV file and a PNG file.
The CSV file contains the raw data, and the PNG file contains a box plot of the differences
between QM and MM values for the specified topology pattern.
"""

import collections
import pathlib
import typing
import click
import tqdm
import json
import time
import sys

from click_option_group import optgroup

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc
import pyarrow.dataset as ds
from openff.toolkit import Molecule

import seaborn as sns
from matplotlib import pyplot as plt

from loguru import logger

logger.remove()
logger.add(sys.stdout)

def batch_filter(
    mapped_smiles: list[str],
    pattern: str = None,
    topology_group: str = None,
) -> dict[str, list[str]]:
    """
    For each SMILES in the input list, find all matches to the given pattern.
    Return a dictionary mapping each SMILES to a list of tuples of atom indices
    that match the pattern.

    Parameters
    ----------
    mapped_smiles : list[str]
        List of mapped SMILES strings to search.
    pattern : str
        The SMARTS pattern to search for.
    topology_group : str
        The topology group to search for, e.g., "Bonds", "Angles", "ImproperTorsions", "ProperTorsions".
    """

    pattern_matches = {}
    for smiles in tqdm.tqdm(mapped_smiles):
        mol = Molecule.from_mapped_smiles(
            smiles,
            allow_undefined_stereo=True
        )
        matches = mol.chemical_environment_matches(pattern)
        mol_matches = []
        if matches:
            for match in matches:
                # assume pattern was ordered sensibly
                if topology_group == "ImproperTorsions":
                    central = match[1]
                    others = sorted([match[0], match[2], match[3]])
                    match = [others[0], central, others[1], others[2]]
                elif match[0] > match[-1]:
                    match = match[::-1]
                mol_matches.append(tuple(match))
            pattern_matches[smiles] = mol_matches
    return pattern_matches

def unwrap_angle(
    reference_angle: float,
    angle: float,
) -> float:
    """Unwrap an angle to be within 180 degrees of a reference angle"""
    while angle - reference_angle > 180:
        angle -= 360
    while angle - reference_angle < -180:
        angle += 360
    return angle

def unwrap_bond(
    reference_bond: float,
    bond: float,
) -> float:
    """Dummy function for bonds"""
    return bond

def compute_bond_difference(mm_value: float, qm_value: float) -> float:
    """Compute difference in bond lengths"""
    return mm_value - qm_value

def compute_angle_difference(mm_value: float, qm_value: float) -> float:
    """
    Compute difference in angles, accounting for periodicity.

    Parameters
    ----------
    mm_value : float
        The angle value from the force field (degrees).
    qm_value : float
        The angle value from QM (degrees).
    
    Returns
    -------
    float
        The difference in angles (degrees), in the range [-180, 180].
    """
    difference = mm_value - qm_value
    if difference > 180:
        difference -= 360
    if difference < -180:
        difference += 360
    return difference


def get_unique_values(
    dataset: ds.Dataset,
    column_name: str = "mapped_smiles",
) -> set[typing.Any]:
    """
    Get the unique values in a column of a dataset in a memory efficient way.

    Parameters
    ----------
    dataset : ds.Dataset
        The dataset to search.
    column_name : str
        The column to get unique values from.

    Returns
    -------
    set
        A set of unique values in the specified column.
    """
    scanner = dataset.scanner(columns=[column_name])
    unique_values = set()
    for batch in scanner.to_batches():
        col = batch.column(column_name)
        unique_values.update(col.to_pylist())
    return unique_values

                    

@click.command(help=__doc__)
@click.option(
    "--input-directory",
    "-i",
    "input_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="topology-values",
    help="Directory containing data files",
)
@click.option(
    "--output-directory",
    "-o",
    "output_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="topology-comparisons",
    help="Directory to write the output files to",
)
@click.option(
    "--input-file",
    "-f",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default="comparison-patterns.json",
    help="JSON file containing the patterns to search for.",
)
@click.option(
    "--index",
    "-n",
    type=int,
    default=0,
    help="Index of the pattern to search for in the input file.",
    show_default=True,
)
@optgroup.group("Parallelization configuration")
@optgroup.option(
    "--n-workers",
    help="The number of workers to distribute the labelling across. Use -1 to request "
    "one worker per batch.",
    type=int,
    default=1,
    show_default=True,
)
@optgroup.option(
    "--worker-type",
    help="The type of worker to distribute the labelling across.",
    type=click.Choice(["lsf", "local", "slurm"]),
    default="local",
    show_default=True,
)
@optgroup.option(
    "--batch-size",
    help="The number of molecules to processes at once on a particular worker.",
    type=int,
    default=500,
    show_default=True,
)
@optgroup.group("Cluster configuration", help="Options to configure cluster workers.")
@optgroup.option(
    "--memory",
    help="The amount of memory (GB) to request per queue worker.",
    type=int,
    default=3,
    show_default=True,
)
@optgroup.option(
    "--walltime",
    help="The maximum wall-clock hours to request per queue worker.",
    type=int,
    default=2,
    show_default=True,
)
@optgroup.option(
    "--queue",
    help="The SLURM queue to submit workers to.",
    type=str,
    default="cpuqueue",
    show_default=True,
)
@optgroup.option(
    "--conda-environment",
    help="The conda environment that SLURM workers should run using.",
    type=str,
)
def main(
    input_file: str,
    index: int,
    input_directory: str = "topology-values",
    output_directory: str = "topology-comparisons",
    worker_type: typing.Literal["slurm", "local"] = "local",
    queue: str = "free",
    conda_environment: str = "ib-dev",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):
    from openff.nagl.utils._parallelization import batch_distributed
    from dask import distributed

    logger.info(f"{time.ctime()} - Starting batch filter")
    start_time = time.time()

    with open(input_file, "r") as f:
        contents = json.load(f)

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)
    image_directory = output_directory / "images"
    image_directory.mkdir(exist_ok=True, parents=True)

    kwargs = contents[index]
    pattern = kwargs["pattern"]
    topology_group = kwargs["topology_group"]
    output_stem = kwargs["output_stem"]

    logger.info(f"Searching for {pattern} and saving to {output_stem}")

    input_directory = pathlib.Path(input_directory)
    input_dataset = ds.dataset(input_directory)
    logger.info(f"Loaded {input_dataset.count_rows()} rows from {input_directory}")

    unique_smiles = get_unique_values(input_dataset, column_name="mapped_smiles")
    logger.info(f"Found {len(unique_smiles)} unique SMILES in {input_directory}")
    unique_smiles = sorted(unique_smiles)

    filtered_smiles_indices = {}

    with batch_distributed(
        unique_smiles,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(
            batch_filter,
            pattern=pattern,
            topology_group=topology_group,
        ))
        for future in tqdm.tqdm(
            distributed.as_completed(futures, raise_errors=False),
            total=len(futures),
            desc="Filtering topology batches",
        ):
            filtered_smiles_indices.update(future.result())

    logger.info(f"Found {len(filtered_smiles_indices)} matches for {pattern} in {topology_group}")

    initial_expression = (
        (pc.field("mapped_smiles").isin(list(filtered_smiles_indices)))
        & (pc.field("topology") == topology_group)
    )


    subset = input_dataset.filter(initial_expression)
    logger.info(f"Filtered to {subset.count_rows()} matching rows in dataset")

    qm_subset = subset.filter(pc.field("method") == "qm")
    qm_rows = qm_subset.to_table(
        columns=["qcarchive_id", "mapped_smiles", "atom_indices", "value"]
    ).to_pylist()
    all_qm_values = collections.defaultdict(dict)
    for qm_row in tqdm.tqdm(qm_rows, desc="Getting QM values"):
        atom_indices = tuple(qm_row["atom_indices"])
        mapped_smiles = qm_row["mapped_smiles"]
        if atom_indices in filtered_smiles_indices[mapped_smiles]:
            all_qm_values[qm_row["qcarchive_id"]][atom_indices] = qm_row["value"]

    ff_subset = subset.filter(pc.field("method") != "qm")
    ff_rows = ff_subset.to_table(
        columns=["qcarchive_id", "mapped_smiles", "atom_indices", "value", "method"]
    ).to_pylist()

    diff_func = compute_bond_difference if topology_group == "Bonds" else compute_angle_difference
    unwrap_func = unwrap_bond if topology_group == "Bonds" else unwrap_angle

    all_output_entries = []
    for ff_row in tqdm.tqdm(ff_rows, desc="Getting FF values"):
        atom_indices = tuple(ff_row["atom_indices"])
        qcarchive_id = ff_row["qcarchive_id"]
        if atom_indices in all_qm_values[qcarchive_id]:
            qm_value = all_qm_values[qcarchive_id][atom_indices]
            unwrapped = unwrap_func(qm_value, ff_row["value"])
            difference = diff_func(ff_row["value"], qm_value)
            check_ = np.isclose(difference, unwrapped - qm_value, atol=1e-4)

            assert check_, f"{difference} != {(unwrapped - qm_value)}"
            entry = {
                "id": f"{qcarchive_id}_{output_stem}_{atom_indices}",
                "qcarchive_id": qcarchive_id,
                "mapped_smiles": ff_row["mapped_smiles"],
                "atom_indices": list(atom_indices),
                "mm_value": ff_row["value"],
                "mm_value_unwrapped": unwrapped,
                "qm_value": qm_value,
                "difference": difference,
                "topology_group": topology_group,
                "method": ff_row["method"].replace("_unconstrained", ""),
            }
            all_output_entries.append(entry)

    df = pd.DataFrame(all_output_entries)
    df["pattern"] = pattern
    df["stem"] = kwargs["output_stem"]
    output_csv = output_directory / "csvs" / (output_stem + ".csv")
    output_csv.parent.mkdir(exist_ok=True, parents=True)
    df.to_csv(output_csv, index=False)
    logger.info(f"Saved {len(df)} to {output_csv}")

    table = pa.Table.from_pandas(df)
    output_pq = output_directory / "pyarrow" / (output_stem + ".parquet")
    output_pq.parent.mkdir(exist_ok=True, parents=True)

    pq.write_table(table, output_pq)
    

    # plot
    ax = sns.boxplot(
        data=df,
        x="difference",
        y="method",
    )
    n = len(df.qcarchive_id.unique())
    ax.set_title(f"{pattern} (n={n})")
    output_png = image_directory / (output_stem + ".png")
    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    logger.info(f"Saved plot to {output_png}")
    

    logger.info(f"{time.ctime()} - Finished getting pattern values")
    elapsed_time = time.time() - start_time
    logger.info(f"Elapsed time: {elapsed_time / 60:.2f} min")
    logger.info("Done!")


if __name__ == "__main__":
    main()

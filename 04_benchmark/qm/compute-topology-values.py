import pathlib
import typing
import click
import tqdm
import time
import sys
from loguru import logger

from click_option_group import optgroup

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

import MDAnalysis as mda
from openff.toolkit import Molecule

logger.remove()
logger.add(sys.stdout)

MDA_TO_OPENFF = {
    "bonds": "Bonds",
    "angles": "Angles",
    "dihedrals": "ProperTorsions",
    "impropers": "ImproperTorsions",
}

def compute_single(
    row,
):
    from openff.units import unit
    mol = Molecule.from_mapped_smiles(
        row["mapped_smiles"],
        allow_undefined_stereo=True
    )
    positions = np.array(row["coordinates"]).reshape((-1, 3))
    mol._conformers = [positions * unit.angstrom]

    base_entry = {
        "qcarchive_id": row["qcarchive_id"],
        "mapped_smiles": row["mapped_smiles"],
        "method": row["method"]
    }

    u = mda.Universe(mol.to_rdkit(), to_guess=["angles", "dihedrals", "impropers"])
    topology_groups = ["bonds", "angles", "dihedrals", "impropers"]

    entries = []
    for topology_group in topology_groups:
        group = getattr(u, topology_group)
        if not len(group):
            continue
        atom_indices = group.indices
        values = group.values()
        if topology_group != "bonds":
            values = np.rad2deg(values)
        for ix, val in zip(atom_indices, values):
            if topology_group == "impropers":
                # reorder the indices to match OpenFF
                central = ix[0]
                others = sorted(ix[1:])
                ix = [others[0], central, others[1], others[2]]

            entry = dict(base_entry)
            entry.update(
                {
                    "topology": MDA_TO_OPENFF[topology_group],
                    "atom_indices": list(ix),
                    "value": val,
                }
            )
            entries.append(entry)
    return entries




def batch_optimize(
    qcarchive_ids: list[str],
    data_directory: str,
):
    dataset = ds.dataset(data_directory)
    subset = dataset.filter(
        pc.field("qcarchive_id").isin(qcarchive_ids)
    )
    rows = subset.to_table(
        columns=["qcarchive_id", "mapped_smiles", "coordinates", "method"]
    ).to_pylist()

    entries = []
    for row in tqdm.tqdm(rows):
        entries.extend(compute_single(row))
    return entries



@click.command()
@click.option(
    "--input-directory",
    "-i",
    "input_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="data",
    help="Directory containing data files",
)
@click.option(
    "--output-directory",
    "-o",
    "output_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="topology-values",
    help="Directory to write output files to",
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
    method_name: str = None,
    input_directory: str = "data",
    output_directory: str = "topology-values",
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

    logger.info(f"{time.ctime()} - Starting batch optimization")
    start_time = time.time()

    input_directory = pathlib.Path(input_directory)
    input_dataset = ds.dataset(input_directory)
    if method_name:
        input_dataset = input_dataset.filter(
            pc.field("method") == method_name
        )
    logger.info(f"Loaded {input_dataset.count_rows()} rows from {input_directory}")
    
    input_qcarchive_ids = input_dataset.to_table(
        columns=["qcarchive_id"]
    ).to_pydict()["qcarchive_id"]
    input_qcarchive_ids = set(input_qcarchive_ids)

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    output_dataset = ds.dataset(output_directory)
    n_files = 0
    if method_name and output_dataset.count_rows():
        output_dataset = output_dataset.filter(
            pc.field("method") == method_name
        )
    if output_dataset.count_rows():
        existing_qcarchive_ids = output_dataset.to_table(
            columns=["qcarchive_id"]
        ).to_pydict()["qcarchive_id"]
        logger.info(f"Loaded {len(existing_qcarchive_ids)} rows from {output_directory}")
        input_qcarchive_ids -= set(existing_qcarchive_ids)
        logger.info(f"Filtered to {len(input_qcarchive_ids)} new rows to process")
        n_files  = len(output_dataset.files)

    input_qcarchive_ids = sorted(input_qcarchive_ids)

    with batch_distributed(
        input_qcarchive_ids,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(
            batch_optimize,
            data_directory=str(input_directory.resolve()),
        ))
        for future in tqdm.tqdm(
            distributed.as_completed(futures, raise_errors=False),
            total=len(futures),
            desc="Calculating topology batches",
        ):
            entries = future.result()
            if len(entries):
                table = pa.Table.from_pylist(entries)
                table_file = output_directory / f"batch-{n_files:04d}.parquet"
                pq.write_table(table, table_file)
                logger.info(f"Wrote {len(entries)} entries to {table_file}")
                n_files += 1

    logger.info(f"{time.ctime()} - Finished batch optimization")
    elapsed_time = time.time() - start_time
    logger.info(f"Elapsed time: {elapsed_time / 60:.2f} min")
    logger.info("Done!")


if __name__ == "__main__":
    main()

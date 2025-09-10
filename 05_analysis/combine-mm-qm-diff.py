"""
Combine multiple MM v. QM difference datasets into a single CSV file.
This also adds labels for which parameters are being labelled.

Note: topology values labelled with parameters do not necessarily match
that parameter in the normal hierarchical force field sense.
Instead, we used the SMARTS to do a SMARTs search.

However, we do try to deduplicate records.
So even if the same bond matches multiple parameter SMARTS,
it is only marked as True in the Bonds column once.
"""
from collections import defaultdict
import pathlib
import sys

import click

from loguru import logger

import numpy as np
import pyarrow.compute as pc
import pyarrow.dataset as ds

logger.remove()
logger.add(sys.stdout)


@click.command(help=__doc__)
@click.option(
    "--input-directory",
    "-i",
    "input_directory",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, readable=True),
    help="Input directory containing a PyArrow dataset with MM v. QM differences.",
    required=True,
)
@click.option(
    "--output-file",
    "-o",
    "output_file",
    type=click.Path(exists=False, dir_okay=False, file_okay=True, writable=True),
    default="output/mm-qm-differences-labelled.csv",
    help="Output CSV file.",
    show_default=True,
)
@click.option(
    "--force-field-col",
    "-cff",
    type=str,
    default="FF",
    help="Column name for force field in the output CSV file.",
)
@click.option(
    "--name-and-forcefield",
    "-nf",
    type=(str, str),
    multiple=True,
    help=(
        "Name and force field pairs to use for naming. "
        "These should be in the format -nf <name> <forcefield>"
    )
)
def main(
    input_directory: str,
    output_file: str,
    force_field_col: str = "FF",
    name_and_forcefield: list[tuple[str, str]] | None = None,
):
    STEM_TO_NAME = {
        v: k for k, v in name_and_forcefield
    }

    dataset = ds.dataset(input_directory)
    logger.info(f"Loaded {dataset.count_rows()} records from {input_directory}")

    # filter to ffs
    dataset = dataset.filter(
        pc.field("method").isin(STEM_TO_NAME.keys())
    )
    logger.info(f"Filtered to {dataset.count_rows()} records with specified force fields")

    df = dataset.to_table().to_pandas()

    # filter only for ids that are present in all methods
    n_ffs = df.method.unique()
    counts = df["id"].value_counts()
    valid_ids = counts[counts == len(n_ffs)].index
    df = df[df["id"].isin(valid_ids)]
    logger.info(f"Filtered to {len(df)} records with all {len(n_ffs)} force fields present")

    # map force field names
    df[force_field_col] = [STEM_TO_NAME.get(v, v) for v in df["method"]]

    # set up label IDs
    unique_parameter_ids = set()
    for stem in df.stem.values:
        _, parameter_id = stem.split("-", 1)
        unique_parameter_ids.add(parameter_id)
    logger.info(f"Found {len(unique_parameter_ids)} unique parameter IDs")

    group_labels = {}
    parameter_types = ["Bonds", "Angles", "ProperTorsions", "ImproperTorsions"]
    for col in parameter_types + sorted(unique_parameter_ids):
        group_labels[col] = np.zeros(len(df), dtype=bool)

    # try to deduplicate to avoid double-counting bonds etc.
    SEEN_TOPOLOGY_GROUPS = defaultdict(set)
    for i, (_, row) in enumerate(df.iterrows()):
        topology_group, parameter_id = row["stem"].split("-", 1)
        key = (row["qcarchive_id"], topology_group, row["atom_indices"])
        if key not in SEEN_TOPOLOGY_GROUPS[topology_group]:
            SEEN_TOPOLOGY_GROUPS[topology_group].add(key)
            group_labels[topology_group][i] = True
        if key not in SEEN_TOPOLOGY_GROUPS[parameter_id]:
            SEEN_TOPOLOGY_GROUPS[parameter_id].add(key)
            group_labels[parameter_id][i] = True

    for col, values in group_labels.items():
        df[col] = values

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_file, index=False)
    logger.info(f"Wrote {len(df)} records to {output_file}")


if __name__ == "__main__":
    main()

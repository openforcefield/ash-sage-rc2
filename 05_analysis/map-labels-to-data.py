"""
This script labels groups in a dataset and identifies groups
overrepresented in a high-error set, selected on a specified metric.
Groups are selected using Checkmol and a force field.

This saves an output CSV that includes all columns of the input data, as well as
labels for the groups identified in each property.
"""


from collections import defaultdict
import pathlib
import sys

import click
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
from rdkit import Chem

from openff.toolkit import Molecule, ForceField

from loguru import logger

logger.remove()
logger.add(sys.stdout)

def sanitize_smiles(smiles: str) -> str:
    """Sanitize SMILES"""
    return Chem.MolToSmiles(
        Chem.MolFromSmiles(smiles)
    )


@click.command(help=__doc__)
@click.option(
    "--input",
    "-i",
    "input_paths",
    multiple=True,
    type=click.Path(exists=True, dir_okay=True, file_okay=True, readable=True),
    help=(
        "Input paths. These can be files or directories. "
        "If a file, it is assumed to be a Pandas CSV. "
        "If a directory, it is assumed to contain a PyArrow dataset. "
        "These are loaded and concatenated into a single dataframe."
    )
)
@click.option(
    "--output",
    "-o",
    "output_file",
    default="output",
    type=click.Path(dir_okay=False, file_okay=True, writable=True),
)
@click.option(
    "--smiles-col",
    "-s",
    "smiles_columns",
    multiple=True,
    help=(
        "SMILES columns to use for grouping. "
        "These should be the names of the columns in the input dataframe."
    )
)
@click.option(
    "--method-col",
    "-cm",
    type=str,
    default="forcefield",
    help="Name of the column containing the method name (usually full name of FF)."
)
@click.option(
    "--force-field-col",
    "-cff",
    type=str,
    default="FF",
    help=(
        "Name of the column containing the human-friendly force field name, "
        "created with `name-and-forcefield`."
    )
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
@click.option(
    "--checkmol-labels-path",
    "-lcm",
    type=str,
    default="labels/checkmol",
    help="Path to the CheckMol labels dataset."
)
@click.option(
    "--forcefield-path",
    "-ff",
    type=str,
    default="../04_benchmark/forcefields/fb-fit-v3-single-mean-k100.offxml",
    help="Path to the force field file."
)
@click.option(
    "--forcefield-labels-path",
    "-lff",
    type=str,
    default="labels/forcefields/v3-k100",
    help="Path to the force field labels dataset."
)
def main(
    input_paths: list[str],
    output_file: str = "labelled-properties.csv",
    smiles_columns: str = [],
    method_col: str = "forcefield",
    force_field_col: str = "FF",
    name_and_forcefield: list[tuple[str, str]] | None = None,
    checkmol_labels_path: str = "labels/checkmol",
    forcefield_path: str = "../04_benchmark/forcefields/fb-fit-v1-single-mean-k100.offxml",
    forcefield_labels_path: str = "labels/forcefields/v1-k100",
):
    dfs = []
    for input_path in input_paths:
        input_path = pathlib.Path(input_path)
        if input_path.is_dir():
            df_ = ds.dataset(input_path).to_table().to_pandas()
        else:
            df_ = pd.read_csv(input_path)
        logger.info(f"Read {len(df_)} entries from {input_path}")
        dfs.append(df_)
    df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Read {len(df)} entries from {len(input_paths)} files")

    if not name_and_forcefield:
        name_and_forcefield = []

    # overwrite with human readable names
    STEM_TO_NAME = {
        v: k for k, v in name_and_forcefield
    }
    df[force_field_col] = [STEM_TO_NAME.get(v, v) for v in df[method_col]]

    # get all smiles to check for
    all_smiles = set()
    for col in smiles_columns:
        all_smiles |= set(df[col].dropna().unique())
    logger.info(f"Found {len(all_smiles)} unique SMILES to label")

    # load labels
    expression = pc.field("smiles").isin(all_smiles)
    checkmol_labels = ds.dataset(checkmol_labels_path).filter(expression)
    checkmol_labels_df = checkmol_labels.to_table().to_pandas()
    logger.info(f"Loaded {len(checkmol_labels_df)} checkmol labels from {checkmol_labels_path}")
    assert all_smiles.issubset(set(checkmol_labels_df["smiles"].unique()))

    forcefield_labels = ds.dataset(forcefield_labels_path).filter(expression)
    forcefield_labels_df = forcefield_labels.to_table().to_pandas()
    logger.info(f"Loaded {len(forcefield_labels_df)} forcefield labels from {forcefield_labels_path}")
    assert all_smiles.issubset(set(forcefield_labels_df["smiles"].unique()))

    # assign to smiles
    CHECKMOL_GROUP_TO_SMILES = defaultdict(set)
    FORCEFIELD_GROUPS_TO_SMILES = defaultdict(set)
    ALL_CHECKMOL_GROUPS = set()
    ALL_FORCEFIELD_GROUPS = defaultdict(set)
    PARAMETER_IDS = {}
    PARAMETER_TYPES = ["Bonds", "Angles", "ProperTorsions", "ImproperTorsions", "vdW"]

    forcefield = ForceField(forcefield_path)

    for smiles, subdf in checkmol_labels_df.groupby("smiles"):
        groups = set(subdf["group"].unique())
        # SMILES_TO_GROUPS_CHECKMOL[smiles] = sorted(groups)
        ALL_CHECKMOL_GROUPS |= groups

        for group in groups:
            CHECKMOL_GROUP_TO_SMILES[group].add(smiles)

    for smiles, subdf in forcefield_labels_df.groupby("smiles"):
        for parameter_type, subsubdf in subdf.groupby("parameter_type"):
            parameter_handler = forcefield.get_parameter_handler(parameter_type)
            parameter_ids = [p.id for p in parameter_handler.parameters]
            PARAMETER_IDS[parameter_type] = parameter_ids

            groups = set(subsubdf["parameter_id"].unique())
            # SMILES_TO_GROUPS_FORCEFIELD[smiles][parameter_type] = sorted(groups)
            ALL_FORCEFIELD_GROUPS[parameter_type] |= groups

            for group in groups:
                FORCEFIELD_GROUPS_TO_SMILES[group].add(smiles)

    ALL_CHECKMOL_GROUPS -= {
        "Hydroxyl" # largely overlaps with Alcohol
    }
    ALL_CHECKMOL_GROUPS = sorted(ALL_CHECKMOL_GROUPS)
    logger.info(f"Found {len(ALL_CHECKMOL_GROUPS)} unique Checkmol groups")

    ALL_FORCEFIELD_GROUPS = {
        ptype: sorted(pids, key=lambda x: PARAMETER_IDS[ptype].index(x))
        for ptype, pids in ALL_FORCEFIELD_GROUPS.items()
    }
    for parameter_type, parameter_ids in ALL_FORCEFIELD_GROUPS.items():
        logger.info(f"Found {len(parameter_ids)} unique {parameter_type} forcefield parameters")

    # set column to presence of group in property
    for col in ALL_CHECKMOL_GROUPS:
        df[col] = df.apply(
            lambda row: any(
                row[smiles_col] in CHECKMOL_GROUP_TO_SMILES[col]
                for smiles_col in smiles_columns
            ),
            axis=1
        )
    df = df.copy() # stop fragmentation warnings
    for parameter_type in PARAMETER_TYPES:
        for col in ALL_FORCEFIELD_GROUPS[parameter_type]:
            df[col] = df.apply(
                lambda row: any(
                    row[smiles_col] in FORCEFIELD_GROUPS_TO_SMILES[col]
                    for smiles_col in smiles_columns
                ),
            axis=1
            )
        df = df.copy() # stop fragmentation warnings

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_file, index=False)
    logger.info(f"Wrote {len(df)} entries to {output_file}")


if __name__ == "__main__":
    main()

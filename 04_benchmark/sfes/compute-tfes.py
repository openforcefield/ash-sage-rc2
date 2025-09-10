"""
This script computes transfer free energies (TFEs)
from aqueous to non-aqueous solvents.
\b
This script is expected to operate on post-processed
CSV files from the `combine-results.py` script, although
it will operate on any CSV files with the expected columns:
    - Id: str
    - Temperature (K): float
    - Pressure (kPa): float
    - Solvent: str
    - Solute: str
    - Value (kcal / mol): float
    - Uncertainty (kcal / mol): float
    - Reference Value (kcal / mol): float
    - Reference Uncertainty (kcal / mol): float
    - Method: str
    - Dataset: str
\b
The output is a single CSV file with the same columns,
with the following changes:
    - Solvent is replaced with Solvent 1 and Solvent 2 columns
    - Solvent 1 is always "O" (water)
    - Solvent 2 is the non-aqueous solvent
    - The Dataset column is set to "TFEs"
    - The Id column is updated to reflect the HFE and SFE IDs

The output file is saved to the specified output path.
"""
import pathlib
import sys
import click
import tqdm

import pandas as pd

from loguru import logger

logger.remove()
logger.add(sys.stdout)


@click.command(help=__doc__)
@click.option(
    "--input-files",
    "-i",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    multiple=True,
    required=True,
    help=(
        "Paths to the input CSV files. "
        "These can be post-processed files from `combine-results.py`. "
        "Multiple files can be specified. "
        "They must have the colummns specified in the docstring."
    )
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(exists=False, dir_okay=False, writable=True),
    default="output/tfes.csv",
    help="Path to the output CSV file.",
    show_default=True,
)
def main(
    input_files: str = [],
    output_file: str = "output/tfes.csv"
):
    dfs = []
    for input_file in input_files:
        df_ = pd.read_csv(input_file)
        dfs.append(df_)
        logger.info(f"Read {len(df_)} rows from {input_file}")

    df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Combined to {len(df)} rows")

    hfes = df[df.Solvent == "O"]
    sfes = df[df.Solvent != "O"]

    hfe_solutes = hfes.Solute.unique()
    logger.info(f"Found {len(hfe_solutes)} unique solutes in HFE data")
    sfe_solutes = sfes.Solute.unique()
    logger.info(f"Found {len(sfe_solutes)} unique solutes in SFE data")
    overlapping_solutes = set(hfe_solutes).intersection(set(sfe_solutes))
    logger.info(f"Found {len(overlapping_solutes)} overlapping solutes")

    tfe_rows = []
    for solute in tqdm.tqdm(overlapping_solutes, desc="Computing TFEs"):
        sfes_for_solute = sfes[sfes.Solute == solute]
        hfes_for_solute = hfes[hfes.Solute == solute]
        method_to_row = {
            row.Method: row for _, row in hfes_for_solute.iterrows()
        }
        # there should be one water HFE per solute
        if not len(method_to_row) == len(hfes_for_solute):
            raise ValueError(
                f"Found {len(hfes_for_solute)} HFE rows for {solute} "
                f"but {len(method_to_row)} unique methods"
            )
        
        for method, method_df in sfes_for_solute.groupby(by="Method"):
            if method not in method_to_row:
                logger.warning(
                    f"Skipping {solute} {method} as no HFE data found"
                )
                continue
            hfe_row = method_to_row[method]
            hfe_ref_value = hfe_row["Reference Value (kcal / mol)"]
            hfe_ref_error = hfe_row["Reference Uncertainty (kcal / mol)"]
            hfe_value = hfe_row["Value (kcal / mol)"]
            hfe_error = hfe_row["Uncertainty (kcal / mol)"]

            for _, sfe_row in method_df.iterrows():
                assert hfe_row["Temperature (K)"] == sfe_row["Temperature (K)"]
                assert hfe_row["Pressure (kPa)"] == sfe_row["Pressure (kPa)"]

                # compute aq -> nonaq
                ref_value = sfe_row["Reference Value (kcal / mol)"] - hfe_ref_value
                value = sfe_row["Value (kcal / mol)"] - hfe_value

                # get errors
                ref_error = (
                    sfe_row["Reference Uncertainty (kcal / mol)"] ** 2 +
                    hfe_ref_error ** 2
                ) ** 0.5
                error = (sfe_row["Uncertainty (kcal / mol)"] ** 2 + hfe_error ** 2) ** 0.5

                row = {
                    "Id": f"{hfe_row.Id} -> {sfe_row.Id}",
                    "Temperature (K)": sfe_row["Temperature (K)"],
                    "Pressure (kPa)": sfe_row["Pressure (kPa)"],
                    "Solvent 1": "O",
                    "Solvent 2": sfe_row["Solvent"],
                    "Solute": solute,
                    "Value (kcal / mol)": value,
                    "Uncertainty (kcal / mol)": error,
                    "Method": method,
                    "Dataset": "TFEs",
                    "Reference Value (kcal / mol)": ref_value,
                    "Reference Uncertainty (kcal / mol)": ref_error,
                }
                tfe_rows.append(row)
        
    df = pd.DataFrame(tfe_rows)
    logger.info(f"Computed {len(df)} TFEs")

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_file, index=False)
    logger.info(f"Wrote TFEs to {output_file}")


if __name__ == "__main__":
    main()

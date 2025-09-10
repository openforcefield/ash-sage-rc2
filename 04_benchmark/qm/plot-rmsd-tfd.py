"""
Plot RMSD and TFD data from QM benchmark.

This script saves six plots:
- rmsd.png: A plot of heavy atom RMSD values for each force field, limited to 1.0 Å.
- rmsd-close.png: A close-up of the RMSD plot, limited to 0.4 Å.
- rmsd-box.png: A box plot of heavy atom RMSD values for each force field.
- tfd.png: A plot of TFD values for each force field, limited to 0.4.
- tfd-close.png: A close-up of the TFD plot, limited to 0.1.
- tfd-box.png: A box plot of TFD values for each force field.
"""

import functools
import pathlib
import logging
import click
import time


import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds


from openff.toolkit import Molecule


import seaborn as sns
from matplotlib import pyplot as plt

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


@functools.cache
def smiles_to_inchi(smiles: str) -> str:
    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    return mol.to_inchi(fixed_hydrogens=True)


def plot_ecdf(
    data: pd.DataFrame,
    x: str,
    palette: dict,
    xmax=1,
    xlabel: str = r"Heavy atom RMSD ($\AA$)",
    imgfile: pathlib.Path = None,
):
    fig, ax = plt.subplots(figsize=(8, 6))
    rmsd_ax = sns.ecdfplot(
        ax=ax,
        data=data,
        x=x,
        hue="Force field",
        palette=palette,
        stat="count",
    )
    rmsd_ax.set_xlim(0, xmax)
    rmsd_ax.set_xlabel(xlabel)
    plt.tight_layout()
    plt.savefig(imgfile, dpi=300)
    logger.info(f"Saved RMSD plot to {imgfile}")
    plt.close()


@click.command()
@click.option(
    "--ff-name-and-stem",
    "-ff",
    "ff_name_and_stems",
    type=(str, str),
    multiple=True,
    help=(
        "Force fields to include and their stem names. "
        "The first argument should be the name of the force field in the plot, "
        "the second argument should be the stem name of the force field (e.g., 'openff_unconstrained-2.2.1'). "
        "This option can be specified multiple times to include multiple force fields."
    )
)
@click.option(
    "--input-directory",
    "-i",
    "input_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="rmsd",
    help="Directory to read RMSD/TFD data from.",
)
@click.option(
    "--output-directory",
    "-o",
    "output_directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="images",
    help="Directory to write output files to.",
)
@click.option(
    "--suffix",
    "-s",
    "suffix",
    type=str,
    default="",
    help="Suffix to append to output file names.",
)
@click.option(
    "--qca-id-col",
    "qcarchive_id_col",
    type=str,
    default="qcarchive_id",
    help="Column to use as QCArchive ID index",
)
def main(
    ff_name_and_stems: list[tuple[str, str]],
    input_directory: str = "rmsd",
    output_directory: str = "images",
    suffix: str = "",
    qcarchive_id_col: str = "qcarchive_id",
):

    logger.info(f"{time.ctime()} - Starting metrics")
    start_time = time.time()

    # load input RMSD/TFD data
    rmsd_dataset = ds.dataset(input_directory)
    FF_STEM_TO_NAME = {
        stem: name for name, stem in ff_name_and_stems
    }
    rmsd_dataset = rmsd_dataset.filter(
        pc.field("method").isin(FF_STEM_TO_NAME.keys())
    )
    rmsd_df = rmsd_dataset.to_table(
        columns=[qcarchive_id_col, "rmsd", "tfd", "method"]
    ).to_pandas()
    logger.info(f"Loaded {len(rmsd_df)} RMSD/TFD records from {input_directory}")

    # count loaded FFs
    unique_ffs = rmsd_df["method"].unique()
    n_unique_ffs = len(unique_ffs)
    logger.info(f"Found {n_unique_ffs} unique force fields in input {input_directory}")

    # filter to only include data points present in all force fields
    counts = rmsd_df.groupby(qcarchive_id_col)["method"].nunique()
    qcarchive_ids = counts[counts == n_unique_ffs].index
    logger.info(f"Found {len(qcarchive_ids)} QCArchive IDs with data in all force fields")
    rmsd_df = rmsd_df[rmsd_df[qcarchive_id_col].isin(qcarchive_ids)]
    rmsd_df["Force field"] = rmsd_df["method"].map(FF_STEM_TO_NAME)

    # set up color palette for plotting
    palette = sns.color_palette(n_colors=len(FF_STEM_TO_NAME))
    ff_names = list(FF_STEM_TO_NAME.values())
    PALETTE = {
        ff_names[i]: palette[i] for i in range(len(FF_STEM_TO_NAME))
    }

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    # plot rmsd, close up, and boxplot
    for xmax, filename in [(1, f"rmsd{suffix}.png"), (0.4, f"rmsd-close{suffix}.png")]:
        plot_ecdf(
            data=rmsd_df,
            x="rmsd",
            palette=PALETTE,
            xmax=xmax,
            xlabel=r"Heavy atom RMSD ($\AA$)",
            imgfile=output_directory / filename,
        )
    ax = sns.boxplot(data=rmsd_df, x="method", hue="Force field", palette=PALETTE, y="rmsd")
    imgfile = output_directory / f"rmsd-box{suffix}.png"
    plt.savefig(imgfile, dpi=300)
    logger.info(f"Saved RMSD box plot to {imgfile}")
    plt.close()

    # plot tfd, close-up, and boxplot

    for xmax, filename in [(0.4, f"tfd{suffix}.png"), (0.1, f"tfd-close{suffix}.png")]:
        plot_ecdf(
            data=rmsd_df,
            x="tfd",
            palette=PALETTE,
            xmax=xmax,
            xlabel="TFD",
            imgfile=output_directory / filename,
        )
    ax = sns.boxplot(data=rmsd_df, x="method", hue="Force field", palette=PALETTE, y="tfd")
    imgfile = output_directory / f"tfd-box{suffix}.png"
    plt.savefig(imgfile, dpi=300)
    logger.info(f"Saved TFD box plot to {imgfile}")
    plt.close()
    logger.info(f"Finished metrics in {time.time() - start_time:.2f} seconds")



if __name__ == "__main__":
    main()

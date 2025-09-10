"""
This script plots the relationship between base and lowest RMSD/TFD values for different force fields as a 2D KDE.
"""
import pathlib

import pyarrow.dataset as ds
import pyarrow.compute as pc
import pandas as pd

from loguru import logger
import click
from matplotlib import pyplot as plt
import seaborn as sns


@click.command(help=__doc__)
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
    "--single-rmsd-dir",
    "-s",
    "single_rmsd_directory",
    default="rmsd-tfd",
    help="Directory containing single/base RMSD/TFD data"
)
@click.option(
    "--all-to-all-rmsd-dir",
    "-a",
    "all_to_all_rmsd_directory",
    default="all-to-all-rmsd-tfd",
    help="Directory containing all-to-all RMSD/TFD data"
)
@click.option(
    "--image-dir",
    "-o",
    "image_directory",
    default="images",
    help="Directory to save output images"
)
def main(
    ff_name_and_stems: list[tuple[str, str]],
    single_rmsd_directory: str = "rmsd-tfd",
    all_to_all_rmsd_directory: str = "all-to-all-rmsd-tfd",
    image_directory: str = "images",
):
    FF_STEM_TO_NAME = {
        stem: name for name, stem in ff_name_and_stems
    }

    # load single QM -> MM data
    single_rmsd_dataset = ds.dataset(single_rmsd_directory)
    logger.info(f"Loaded {single_rmsd_dataset.count_rows()} single RMSD/TFD records from {single_rmsd_directory}")
    single_rmsd_dataset = single_rmsd_dataset.filter(
        pc.field("method").isin(FF_STEM_TO_NAME.keys())
    )
    logger.info(f"Filtered to {single_rmsd_dataset.count_rows()} records with specified force fields")

    # load all-to-all matching where each MM conf is matched to nearest QM conf
    all_to_all_rmsd_dataset = ds.dataset(all_to_all_rmsd_directory)
    logger.info(f"Loaded {all_to_all_rmsd_dataset.count_rows()} all-to-all RMSD/TFD records from {all_to_all_rmsd_directory}")
    all_to_all_rmsd_dataset = all_to_all_rmsd_dataset.filter(
        pc.field("method").isin(FF_STEM_TO_NAME.keys())
    )
    logger.info(f"Filtered to {all_to_all_rmsd_dataset.count_rows()} records with specified force fields")

    single_rmsd_df = single_rmsd_dataset.to_table().to_pandas()
    all_to_all_rmsd_df = all_to_all_rmsd_dataset.to_table().to_pandas()

    assert len(single_rmsd_df) == len(all_to_all_rmsd_df), "Single and all-to-all RMSD DataFrames must have the same length"
    
    # Merge on method + qcarchive_id (from x) == ff_qcarchive_id (from y)
    merged_rmsd = pd.merge(
        all_to_all_rmsd_df,
        single_rmsd_df,
        left_on=["method", "ff_qcarchive_id"],
        right_on=["method", "qcarchive_id"],
        how="inner",
        suffixes=("_lowest", "_raw"),
    )

    # filter on QCArchive IDs being present in all methods
    n_expected_ffs = len(FF_STEM_TO_NAME)

    n_unique_qca_ids = merged_rmsd["qcarchive_id"].nunique()
    counts = merged_rmsd.groupby("qcarchive_id").size()
    valid_ids = counts[counts == n_expected_ffs].index
    n_valid_qca_ids = valid_ids.nunique()
    logger.info(f"Found {n_valid_qca_ids} QCArchive IDs present in all FFs out of {n_unique_qca_ids}")
    merged_rmsd = merged_rmsd[merged_rmsd["qcarchive_id"].isin(valid_ids)]
    logger.info(f"Filtered to {len(merged_rmsd)} observations")

    merged_rmsd["Force field"] = merged_rmsd["method"].map(FF_STEM_TO_NAME)

    # plot
    image_directory = pathlib.Path(image_directory)
    image_directory.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax = sns.kdeplot(
        data=merged_rmsd,
        x="rmsd_lowest",
        y="rmsd_raw",
        hue="Force field",
    )
    ax.set_xlabel(r"Heavy-atom RMSD (Nearest) [$\AA$]")
    ax.set_ylabel(r"Heavy-atom RMSD (Base) [$\AA$]")
    lim = (-0.1, 1.5)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.plot(lim, lim, ls="--", lw=1, color="k")
    plt.tight_layout()
    outfile = image_directory / "rmsd-lowest-to-raw.png"
    plt.savefig(outfile, dpi=300)
    logger.info(f"Saved plot to {outfile}")
    plt.close()

    fig, ax = plt.subplots(figsize=(5, 5))
    ax = sns.kdeplot(
        data=merged_rmsd,
        x="tfd_lowest",
        y="tfd_raw",
        hue="Force field",
    )
    ax.set_xlabel("TFD (Nearest)")
    ax.set_ylabel("TFD (Base)")
    lim = (-0.05, 0.2)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.plot(lim, lim, ls="--", lw=1, color="k")
    plt.tight_layout()
    outfile = image_directory / "tfd-lowest-to-raw.png"
    plt.savefig(outfile, dpi=300)
    logger.info(f"Saved plot to {outfile}")
    plt.close()


if __name__ == "__main__":
    main()

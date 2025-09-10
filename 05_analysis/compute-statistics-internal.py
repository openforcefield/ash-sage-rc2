"""
Calculate statistics for internal coordinates.
This has to be structured a bit differently to the plain `compute-statistics.py`
for size reasons.
\b
This script saves a summary CSV file with the computed statistics.
The CSV contains the columns:
    - force field column (usually 'FF') (str): force field
    - stat (str): statistic in question (RMSE, MUE, MSE, rho, or R2)
    - n (int): the number of observations
    - mle (float): the MLE value
    - mean (float): the mean value
    - low (float): the lower bound of the confidence interval
    - high (float): the upper bound of the confidence interval
    - group (str): the functional group
"""

import click
import pathlib
import sys
import typing

import tqdm
import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds

from loguru import logger
import sklearn.metrics
import scipy.stats

logger.remove()
logger.add(sys.stdout)

def bootstrap_statistic(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    dy_true: np.ndarray | None = None,
    dy_pred: np.ndarray | None = None,
    ci: float = 0.95,
    statistic: str = "RMSE",
    nbootstrap: int = 1000,
    include_true_uncertainty: bool = False,
    include_pred_uncertainty: bool = False,
) -> dict:
    """Compute mean and confidence intervals of specified statistic.

    Copied from cinnabar.stats.bootstrap_statistic, with additions.

    Parameters
    ----------
    y_true : ndarray with shape (N,)
        True values
    y_pred : ndarray with shape (N,)
        Predicted values
    dy_true : ndarray with shape (N,) or None
        Errors of true values. If None, the values are assumed to have no errors
    dy_pred : ndarray with shape (N,) or None
        Errors of predicted values. If None, the values are assumed to have no errors
    ci : float, optional, default=0.95
        Interval for confidence interval (CI)
    statistic : str
        Statistic, one of ['RMSE', 'MUE', 'R2', 'rho','KTAU','RAE', 'MSE']
    nbootstrap : int, optional, default=1000
        Number of bootstrap samples
    include_true_uncertainty : bool, default False
        whether to account for the uncertainty in y_true when bootstrapping
    include_pred_uncertainty : bool, default False
        whether to account for the uncertainty in y_pred when bootstrapping

    Returns
    -------
    rmse_stats : dict of float
        'mean' : mean RMSE
        'stderr' : standard error
        'low' : low end of CI
        'high' : high end of CI
    """

    def compute_statistic(
        y_true_sample: np.ndarray,
        y_pred_sample: np.ndarray,
        statistic: str
    ) -> float:
        """Compute requested statistic.

        Parameters
        ----------
        y_true : ndarray with shape (N,)
            True values
        y_pred : ndarray with shape (N,)
            Predicted values
        statistic : str
            Statistic, one of ['RMSE', 'MUE', 'R2', 'rho','RAE','KTAU']

        """

        def calc_RAE(y_true_sample: np.ndarray, y_pred_sample: np.ndarray):
            MAE = sklearn.metrics.mean_absolute_error(y_true_sample, y_pred_sample)
            mean = np.mean(y_true_sample)
            MAD = np.sum([np.abs(mean - i) for i in y_true_sample]) / float(len(y_true_sample))
            return MAE / MAD

        def calc_RRMSE(y_true_sample: np.ndarray, y_pred_sample: np.ndarray):
            rmse = np.sqrt(sklearn.metrics.mean_squared_error(y_true_sample, y_pred_sample))
            mean_exp = np.mean(y_true_sample)
            mds = np.sum([(mean_exp - i) ** 2 for i in y_true_sample]) / float(len(y_true_sample))
            rrmse = np.sqrt(rmse**2 / mds)
            return rrmse

        if statistic == "RMSE":
            return sklearn.metrics.root_mean_squared_error(y_true_sample, y_pred_sample)
        elif statistic == "MUE":
            return sklearn.metrics.mean_absolute_error(y_true_sample, y_pred_sample)
        # mean signed error
        elif statistic == "MSE":
            return np.mean(y_pred_sample - y_true_sample)
        elif statistic == "R2":
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(y_true_sample, y_pred_sample)
            return r_value**2
        elif statistic == "rho":
            return scipy.stats.pearsonr(y_true_sample, y_pred_sample)[0]
        elif statistic == "RAE":
            return calc_RAE(y_true_sample, y_pred_sample)
        elif statistic == "KTAU":
            return scipy.stats.kendalltau(y_true_sample, y_pred_sample)[0]
        else:
            raise Exception("unknown statistic '{}'".format(statistic))

    assert len(y_true) == len(y_pred)

    sample_size = len(y_true)
    s_n = np.zeros([nbootstrap], np.float64)  # s_n[n] is the statistic computed for bootstrap sample n
    for replicate in tqdm.tqdm(range(nbootstrap)):
        js = np.random.choice(np.arange(sample_size), size=[sample_size], replace=True)
        y_true_sample = np.random.normal(
            loc=y_true[js],
            scale=0
        )
        y_pred_sample = np.random.normal(
            loc=y_pred[js],
            scale=0
        )
        assert len(y_true_sample) == len(y_pred_sample) == sample_size
        # y_true_sample = np.zeros_like(y_true)
        # y_pred_sample = np.zeros_like(y_pred)
        # for i, j in enumerate(np.random.choice(np.arange(sample_size), size=[sample_size], replace=True)):
        #     y_true_sample[i] = np.random.normal(loc=y_true[j], scale=0)
        #     y_pred_sample[i] = np.random.normal(loc=y_pred[j], scale=0)
        try:
            s_n[replicate] = compute_statistic(y_true_sample, y_pred_sample, statistic)
        except Exception as e:
            print(y_true_sample, y_pred_sample)
            raise e

    rmse_stats = dict()
    rmse_stats["mle"] = compute_statistic(y_true, y_pred, statistic)
    rmse_stats["stderr"] = np.std(s_n)
    rmse_stats["mean"] = np.mean(s_n)
    # TODO: Is there a canned method to do this?
    s_n = np.sort(s_n)
    low_frac = (1.0 - ci) / 2.0
    high_frac = 1.0 - low_frac
    rmse_stats["low"] = s_n[int(np.floor(nbootstrap * low_frac))]
    rmse_stats["high"] = s_n[int(np.ceil(nbootstrap * high_frac))]

    return rmse_stats



def compute_stats(
    df: pd.DataFrame,
    forcefield_col: str = "FF",
    x_col: str = None,
    y_col: str = None,
) -> pd.DataFrame:
    stat_rows = []

    for ff, subdf in df.groupby(by=forcefield_col):
        y = subdf[y_col].values
        x = subdf[x_col].values

        for i, stat in enumerate(["RMSE", "MUE", "MSE", "R2", "rho"], 1):
            s = bootstrap_statistic(
                x,
                y,
                statistic=stat,
                include_true_uncertainty=True,
                include_pred_uncertainty=True,
            )
            stat_row = {forcefield_col: ff, "stat": stat, "n": len(subdf)}
            stat_row.update(s)
            stat_rows.append(stat_row)
    return pd.DataFrame(stat_rows)


@click.command(help=__doc__)
@click.option(
    "--input",
    "-i",
    "input_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Input directory containing data calculate stats for.",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    help="Output CSV file for saving computed statistics.",
)
@click.option(
    "--x-col",
    "-cx",
    "x_col",
    type=str,
    default="qm_value",
    help="Column name for x-axis."
)
@click.option(
    "--y-col",
    "-cy",
    "y_col",
    type=str,
    default="mm_value_unwrapped",
    help="Column name for y-axis."
)
@click.option(
    "--topology-group",
    "-g",
    "topology_group",
    type=click.Choice(["Bonds", "Angles", "ProperTorsions", "ImproperTorsions"], case_sensitive=True),
    required=True,
    help="Topological group to compute statistics for.",
)
@click.option(
    "--forcefield-col",
    "-cff",
    "forcefield_col",
    type=str,
    default="FF",
    help="Column name for forcefield.",
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
    x_col: str = None,
    y_col: str = None,
    topology_group: typing.Literal["Bonds", "Angles", "ProperTorsions", "ImproperTorsions"] | None = None,
    forcefield_col: str = "FF",
    name_and_forcefield: list[tuple[str, str]] | None = None,
):
    
    STEM_TO_NAME = {}
    if name_and_forcefield:
        STEM_TO_NAME = {
            v: k for k, v in name_and_forcefield
        }

    dataset = ds.dataset(input_directory).filter(
        pc.field("topology_group") == topology_group
    )
    logger.info(f"Loaded dataset with {dataset.count_rows()} rows")

    # partition on stem
    output_stems = dataset.to_table(columns=["stem"]).to_pydict()["stem"]
    groups_to_compute = sorted(set(output_stems))
    logger.info(f"Found {len(groups_to_compute)} groups to compute")

    # assert topology_group in groups_to_compute
    # # move to start
    # groups_to_compute.remove(topology_group)
    # groups_to_compute.insert(0, topology_group)

    if not STEM_TO_NAME:
        methods = dataset.to_table(columns=["method"]).to_pydict()["method"]
        STEM_TO_NAME = {v: v for v in set(methods)}

    dataset = dataset.filter(
        pc.field("method").isin(STEM_TO_NAME.keys())
    )
    logger.info(f"Filtered to {dataset.count_rows()} records with specified force fields")

    stat_dfs = []
    for group in tqdm.tqdm(groups_to_compute, desc="Computing stats"):
        sub_dataset = dataset.filter(pc.field("stem") == group)
        df = sub_dataset.to_table(
            columns=["id", "method", x_col, y_col]
        ).to_pandas()
        logger.info(f"Group {group}: found {len(df)} entries")

        # restrict observations to those with all force fields present
        n_ffs = df.method.unique()
        counts = df["id"].value_counts()
        valid_ids = counts[counts == len(n_ffs)].index
        df = df[df["id"].isin(valid_ids)]
        logger.info(f"Group {group}: filtered to {len(df)} records with all {len(n_ffs)} force fields present")

        df[forcefield_col] = [STEM_TO_NAME.get(v, v) for v in df["method"]]

        stat_df = compute_stats(
            df, forcefield_col=forcefield_col,
            x_col=x_col, y_col=y_col
        )
        if group == topology_group:
            stat_df["group"] = "All"
        else:
            # just use parameter ID/name
            stat_df["group"] = group.split("-", 1)[-1]
        stat_dfs.append(stat_df)

    # Concatenate all stats DataFrames
    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if not (output_file.is_file() and output_file.exists()):
        all_stats_df = pd.concat(stat_dfs, ignore_index=True)
        all_stats_df.to_csv(output_file, index=False)
        logger.info(f"Saved combined stats to {output_file}")


if __name__ == "__main__":
    main()

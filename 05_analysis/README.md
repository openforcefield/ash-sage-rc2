# Procedure

1. Run `label-with-checkmol.py` and `label-with-forcefield.py` to get labels for unique SMILES in a particular dataset (see `run-label.sh`)
2. Run `map-labels-to-data.sh` to map those labels to each property in a dataset.
3. Run `run-compute-statistics.sh` to compute statistics.
4. Run `run-plot-aggregated-statistics.sh` to plot aggregated statistics
5. Run `plot-labelled-scatter.py` to plot scatter plots instead of stats
6. Run `run-get-largest-differences.sh` to get the largest meaningful differences in a dataset. Results are written to logs
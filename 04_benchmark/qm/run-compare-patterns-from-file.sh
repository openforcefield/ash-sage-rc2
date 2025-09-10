#!/usr/bin/env bash
#SBATCH -J compare-patterns
#SBATCH --array=1-3
#SBATCH -p standard
#SBATCH -t 16:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=196gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A-%a.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate nagl-valence

export OE_LICENSE="/data/homezvol3/lilyw7/oe_license.txt"


echo $SLURM_ARRAY_TASK_ID

python compare-topology-pattern-from-file.py                    \
        --n-workers                     300                     \
        --worker-type                   "slurm"                 \
        --batch-size                    500                     \
        --memory                        8                       \
        --walltime                      480                     \
        --queue                         "free"                  \
        --conda-environment             "nagl-valence"          \
    -i       "topology-values"                                  \
    -f       "comparison-patterns/base-patterns.json"  \
    -n       $SLURM_ARRAY_TASK_ID                               \
    -o       "topology-comparisons"


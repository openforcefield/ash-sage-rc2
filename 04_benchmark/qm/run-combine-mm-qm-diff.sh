#!/usr/bin/env bash
#SBATCH -J combine-qm-mm
#SBATCH -p standard
#SBATCH -t 16:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=164gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate nagl-valence

QM_FFS=(
    -nf 'Sage 2.0.0' 'openff-2.0.0'
    -nf 'Sage 2.2.1' 'openff-2.2.1'
    -nf 'Sage 2.3.0rc1' 'openff-2.3.0rc1'
    -nf 'v1-k100'    'fb-fit-v1-single-mean-k100'
    -nf 'v3-k100'    'fb-fit-v3-single-mean-k100'
)
python combine-mm-qm-diff.py \
    "${QM_FFS[@]}" \
    -i  topology-comparisons/pyarrow \
    -o  output/qm/mm-qm-diff/labelled-properties.csv > logs/combine-mm-qm-diff.log


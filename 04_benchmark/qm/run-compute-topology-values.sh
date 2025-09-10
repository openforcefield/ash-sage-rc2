#!/usr/bin/env bash
#SBATCH -J compute-topology-values
#SBATCH -p standard
#SBATCH -t 16:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate nagl-valence

export OE_LICENSE="/data/homezvol3/lilyw7/oe_license.txt"

FFNAME=qm

# compute topology values for ICRMSD comparison
python compute-topology-values.py                              \
       --n-workers                     300                     \
       --worker-type                   "slurm"                 \
       --batch-size                    500                     \
       --memory                        8                       \
       --walltime                      480                     \
       --queue                         "free"                  \
       --conda-environment             "nagl-valence"          \
   -i       "data/optimization/${FFNAME}"                      \
   -o       "topology-values/${FFNAME}" > logs/compute-topology-values-$FFNAME.log



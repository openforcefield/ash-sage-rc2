#!/usr/bin/env bash
#SBATCH -J get-a2a-rmsd-tfd
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

FFDIR=""
FFNAME="openff_unconstrained-2.2.1"
FFNAME="openff_unconstrained-2.1.0"
FFNAME="openff_unconstrained-2.0.0"
FFNAME="openff_unconstrained-1.3.1"

FFDIR="../forcefields/"
FFNAME=fb-fit-v0-single-mean-k20_unconstrained
FFNAME=fb-fit-v0-single-mean-k100_unconstrained
FFNAME=fb-fit-v1-single-mean-k100_unconstrained
FFNAME=fb-fit-v2-single-mean-k100_unconstrained
FFNAME=fb-fit-v3-single-mean-k100_unconstrained
FFNAME='openff_unconstrained-2.3.0rc1'
FFNAME="openff_unconstrained-2.2.1-ashgc"


FORCEFIELD="${FFDIR}${FFNAME}.offxml"

echo "Benchmarking $FORCEFIELD"



# get all to all RMSD matrix for ddEs
python get-all-to-all-rmsds-and-tfds.py                         \
        --n-workers                     300                     \
        --worker-type                   "slurm"                 \
        --batch-size                    200                     \
        --memory                        8                       \
        --walltime                      480                     \
        --queue                         "free"                  \
        --conda-environment             "nagl-valence"          \
    -d       "data/optimization"                                \
    -r       "all-to-all-rmsd-tfd"                     \
    -ff      $FORCEFIELD > logs/get-all-to-all-${FFNAME}.log



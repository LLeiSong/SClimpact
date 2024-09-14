#!/bin/bash
#SBATCH -c 10
#SBATCH --mem 180G
#SBATCH -p batch
#SBATCH -t 96:00:00
#SBATCH --job-name=sdm_%a
#SBATCH --output=sdm_%a.out
#SBATCH --error=sdm_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lsong@ucsb.edu
#SBATCH --array=1-145

source ~/anaconda3/bin/activate
conda activate sdm
cd ~

unset DISPLAY
srun Rscript scripts/fit_rf.R -i $SLURM_ARRAY_TASK_ID -c 10

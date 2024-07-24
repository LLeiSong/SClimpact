#!/bin/bash
#SBATCH -c 10
#SBATCH --mem 100G
#SBATCH -p batch
#SBATCH -t 10:00:00
#SBATCH --job-name=maxent_%a
#SBATCH --output=maxent_%a.out
#SBATCH --error=maxent_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lsong@ucsb.edu
#SBATCH --array=1-145

source ~/anaconda3/bin/activate
conda activate sdm
cd ~

unset DISPLAY
srun Rscript scripts/fit_maxent.R -i $SLURM_ARRAY_TASK_ID

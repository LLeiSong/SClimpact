#!/bin/bash
#SBATCH -c 4
#SBATCH --mem 20G
#SBATCH -p largemem
#SBATCH -t 10:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong@ucsb.edu

SP=$1

conda activate sdm
cd ~/SCImpact

srun Rscript R/rf_warp.R -s $SP

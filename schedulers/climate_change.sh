#!/bin/bash
#SBATCH -c 10
#SBATCH --mem 600G
#SBATCH -p largemem
#SBATCH -t 10-00:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong@ucsb.edu

FEATURE=$1

conda activate sdm
cd ~/SCImpact

srun Rscript R/cc_warp.R -f $FEATURE

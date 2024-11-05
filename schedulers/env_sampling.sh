#!/bin/bash
#SBATCH -c 40
#SBATCH --mem 1200G
#SBATCH -p largemem
#SBATCH -t 01-00:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong@ucsb.edu

SP=$1

conda activate sdm
cd ~/SCImpact

srun Rscript R/es_warp.R -s $SP
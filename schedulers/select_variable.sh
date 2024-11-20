#!/bin/bash
#SBATCH -c 3
#SBATCH --mem 30G
#SBATCH -p gpu
#SBATCH --exclude=node112,node113
#SBATCH -t 01-00:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong@ucsb.edu

SP=$1

conda activate sdm
cd ~/SCImpact

srun Rscript R/vs_warp.R -s $SP

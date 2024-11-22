#!/bin/bash
#SBATCH -c 40
#SBATCH --mem 180G
#SBATCH -t 07-00:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong@ucsb.edu

SP=$1

conda activate sdm
cd ~/SCImpact

srun Rscript R/shap_warp.R -s $SP

#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --job-name=VDJ_PLL
#SBATCH --error=log/VDJ_pll_combine.error
#SBATCH --output=log/VDJ_pll_combine.out

cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
#datasets=("OVA_V7" "horns2020a__VDJ_RAW")
datasets="Bruhn"
for data in "${datasets[@]}"
  do
      Rscript vdj_pll_combine.R $data
  done
echo "finished"


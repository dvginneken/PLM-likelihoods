#!/bin/bash
#SBATCH --job-name=pll_concatenate
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --time=1:00:00
#SBATCH --output=log/pll_concatenate.out
#SBATCH --error=log/pll_concatenate.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
#datasets=("OVA_V7" "horns2020a__VDJ_RAW")
datasets="Bruhn"
for data in "${datasets[@]}"
  do
      Rscript pll_concatenate.R $data
  done
echo "finished"

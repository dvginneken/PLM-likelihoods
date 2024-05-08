#!/bin/bash
#SBATCH --job-name=pll_correlations
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --output=log/pll_correlations.out
#SBATCH --error=log/pll_correlations.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
#datasets=("OVA_V7" "horns2020a__VDJ_RAW")
datasets="Bruhn"
for data in "${datasets[@]}"
  do
      Rscript SourceCorrelation_chains.R $data
      Rscript PLMCorrelation_chains.R $data
  done
echo "finished"

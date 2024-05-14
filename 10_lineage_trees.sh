#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --job-name=AntibodyForests
#SBATCH --error=log/LineageTrees.error
#SBATCH --output=log/LineageTrees.out

cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn")
for data in "${datasets[@]}"
  do
      Rscript LineageTrees.R $data
      echo "$data done"
  done
echo "finished"


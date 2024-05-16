#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --job-name=MP_AntibodyForests
#SBATCH --error=log/MP_LineageTrees.error
#SBATCH --output=log/MP_LineageTrees.out

cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn" "Bieberich")
for data in "${datasets[@]}"
  do
      Rscript LineageTrees_mp.R $data
      echo "$data done"
  done
echo "finished"


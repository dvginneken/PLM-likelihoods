#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --job-name=VDJ
#SBATCH --error=log/CreateVDJ.error
#SBATCH --output=log/CreateVDJ.out

cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
#datasets=("OVA_V7" "horns2020a__VDJ_RAW")
datasets="Bruhn"
for data in "${datasets[@]}"
  do
      Rscript CreateVDJ.R $data
  done
echo "finished"


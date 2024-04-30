#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=VDJ
#SBATCH --error=log/CreateVDJ.error
#SBATCH --output=log/CreateVDJ.out

cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW")
for data in "${datasets[@]}"
  do
      Rscript CreateVDJ.R $data
  done
echo "finished"


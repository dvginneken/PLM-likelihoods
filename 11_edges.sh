#!/bin/bash
#SBATCH --job-name=edges
#SBATCH --nodes=1
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --output=log/edges.out
#SBATCH --error=log/edges.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bieberich" "Bruhn")
for data in "${datasets[@]}"
  do
      Rscript transform_edge_dataframe.R $data "MP"
      Rscript transform_edge_dataframe.R $data "default"
  done
echo "finished"

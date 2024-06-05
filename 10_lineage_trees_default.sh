#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G
#SBATCH --partition=cpu
#SBATCH --job-name=default_AntibodyForests
#SBATCH --error=log/default_LineageTrees.error
#SBATCH --output=log/default_LineageTrees.out

cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
#datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn" "Bieberich")
datasets="Kim"
for data in "${datasets[@]}"
do
    Rscript LineageTrees_default.R $data
    echo "$data done"
done
echo "finished"


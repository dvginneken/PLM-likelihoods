#!/bin/bash
#SBATCH --job-name=evovelocity_adata
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=100:00:00
#SBATCH --output=log/evovelocity_adata.out
#SBATCH --error=log/evovelocity_adata.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn")
for data in "${datasets[@]}"
do
    python3 evo_velo_adata_maker.py --dataset=$data --group="sample_id" --include_germline
    echo "done $data"
done
echo "finished"

#!/bin/bash
#SBATCH --job-name=evovelocity
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=100:00:00
#SBATCH --output=log/evovelocity_germlines.out
#SBATCH --error=log/evovelocity_germlines.error

export TORCH_HOME=/hpc/dla_lti/dvanginneken/cache/
cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn" "Kim")
#datasets="OVA_V7"
models=("esm" "protbert" "sapiens" "ablang")
for data in "${datasets[@]}"
do
    for model in "${models[@]}"
    do
        python3 evo_velo_adata_maker_germlines.py --dataset=$data --model=$model
    done
    echo "done $data"
done
echo "finished"



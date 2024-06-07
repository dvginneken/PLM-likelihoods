#!/bin/bash
#SBATCH --job-name=germline_emb
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=100:00:00
#SBATCH --output=log/germline_emb_human.out
#SBATCH --error=log/germline_emb_human.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
#datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn")
datasets="Human"
for data in "${datasets[@]}"
do
    python3 embeddings_calculator.py --dataset=$data --mode="general" --file_path="/hpc/dla_lti/dvanginneken/PLM-likelihoods/data/Human_IMGT_germlines.csv"
    echo "done $data"
done
echo "finished"

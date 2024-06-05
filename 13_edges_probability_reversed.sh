#!/bin/bash
#SBATCH --job-name=prob_edges
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=24gb
#SBATCH --time=10:00:00
#SBATCH --output=log/edges_probability_reversed.out
#SBATCH --error=log/edges_probability_reversed.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Kim" "Bruhn")
#datasets="horns2020a__VDJ_RAW"
for data in "${datasets[@]}"
do
    python3 clonal_evolution_prob_reversed.py -d $data -m "MP"
    python3 clonal_evolution_prob_reversed.py -d $data -m "default"
done
echo "finished"

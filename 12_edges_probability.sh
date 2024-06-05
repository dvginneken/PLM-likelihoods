#!/bin/bash
#SBATCH --job-name=prob_edges
#SBATCH --nodes=1
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1
#SBATCH --mem=24gb
#SBATCH --time=10:00:00
#SBATCH --output=log/edges_probability.out
#SBATCH --error=log/edges_probability.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Kim" "Bruhn")
#datasets="Kim"
for data in "${datasets[@]}"
  do
      python3 clonal_evolution_prob.py -d $data -m "MP"
      python3 clonal_evolution_prob.py -d $data -m "default"
  done
echo "finished"

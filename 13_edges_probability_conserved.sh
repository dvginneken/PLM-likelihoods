#!/bin/bash
#SBATCH --job-name=prob_edges_conserved
#SBATCH --nodes=1
#SBATCH --mem=24gb
#SBATCH --time=4:00:00
#SBATCH --output=log/edges_probability_conserved.out
#SBATCH --error=log/edges_probability_conserved.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bieberich" "Bruhn")
for data in "${datasets[@]}"
  do
      python3 clonal_evolution_prob_conserved.py -d $data -m "MP"
      python3 clonal_evolution_prob_conserved.py -d $data -m "default"
  done
echo "finished"

#!/bin/bash
#SBATCH --job-name=evovelocity_plot
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --output=log/evovelocity_plot.out
#SBATCH --error=log/evovelocity_plot.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn")
for data in "${datasets[@]}"
do
    python3 evo_velo_plotting.py --dataset=$data --group="sample_id" --color="SHM_count" --samplewise --include_germline
    echo "done $data"
done
echo "finished"

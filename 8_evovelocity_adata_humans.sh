#!/bin/bash
#SBATCH --job-name=seqHum_evovelocity
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=100:00:00
#SBATCH --output=log/evovelocity_sequences_human.out
#SBATCH --error=log/evovelocity_sequences_human.error

export TORCH_HOME=/hpc/dla_lti/dvanginneken/cache/
cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts

models=("protbert" "sapiens" "ablang")
#models="esm"
for model in "${models[@]}"
do
    python3 evo_velo_adata_maker_human.py --model=$model --input_source="full_VDJ"
done

echo "finished"



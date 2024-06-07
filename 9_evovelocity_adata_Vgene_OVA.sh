#!/bin/bash
#SBATCH --job-name=evovelocity
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=1:00:00
#SBATCH --output=log/evovelocity_vgene_ova.out
#SBATCH --error=log/evovelocity_vgene_ova.error

export TORCH_HOME=/hpc/dla_lti/dvanginneken/cache/
cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
#models=("esm" "protbert" "sapiens" "ablang")
models=("esm")
for model in "${models[@]}"
do
    python3 evo_velo_adata_maker_OVA_subfamily.py --model=$model --v_gene="IGHV3-1"
    python3 evo_velo_adata_maker_OVA_subfamily.py --model=$model --v_gene="IGHV5-6"
    python3 evo_velo_adata_maker_OVA_subfamily.py --model=$model --v_gene="IGHV9-3"
done
echo "finished"



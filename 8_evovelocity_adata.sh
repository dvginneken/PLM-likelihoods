#!/bin/bash
#SBATCH --job-name=evovelocity
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=100:00:00
#SBATCH --output=log/evovelocity.out
#SBATCH --error=log/evovelocity.error

export TORCH_HOME=/hpc/dla_lti/dvanginneken/cache/
cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
#datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn" "Bieberich")
datasets="Bruhn"
for data in "${datasets[@]}"
do
    python3 evo_velo_adata_maker_karlis.py --dataset=$data --model="esm" --sample_id="Bruhn" --input_source="full_VDJ" --color="c_gene"
    echo "done $data"
done
echo "finished"



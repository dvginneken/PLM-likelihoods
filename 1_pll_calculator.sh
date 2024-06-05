#!/bin/bash
#SBATCH --job-name=pll_calculator
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=50:00:00
#SBATCH --output=log/pll_calculator.out
#SBATCH --error=log/pll_calculator.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
modes=("cdr3_only" "full_VDJ" "cdr3_from_VDJ")
#datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn")
datasets="Kim"
for mode in "${modes[@]}"
do
    for data in "${datasets[@]}"
    do
        python3 pll_calculator.py --dataset=$data --mode=$mode
        echo "done $data $mode"
    done
done
echo "finished"

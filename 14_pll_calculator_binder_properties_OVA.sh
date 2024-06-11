#!/bin/bash
#SBATCH --job-name=pll_calculator
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=50:00:00
#SBATCH --output=log/pll_calculator_binder_properties_OVA.out
#SBATCH --error=log/pll_calculator_binder_properties_OVA.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
python3 pll_calculator.py --dataset="OVA_V7" --mode="general" --file_path="/hpc/dla_lti/dvanginneken/PLM-likelihoods/data/OVA_V7/Binder_properties_cleaned.csv"
echo "finished"

#!/bin/bash
#SBATCH --job-name=pll_calculator
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=50:00:00
#SBATCH --output=log/pll_calculator_affinity_kim.out
#SBATCH --error=log/pll_calculator_affinity_kim.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
python3 pll_calculator.py --dataset="Kim" --mode="general" --file_path="/hpc/dla_lti/dvanginneken/PLM-likelihoods/data/Kim/Affinity/Affinity_dataframe_HC.csv"
echo "finished"

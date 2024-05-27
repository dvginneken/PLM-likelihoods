#!/bin/bash
#SBATCH --job-name=pll_calculator
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=50:00:00
#SBATCH --output=log/pll_calculator_kim.out
#SBATCH --error=log/pll_calculator_kim.error


cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/scripts
python3 pll_calculator.py --dataset="Kim" --mode="general" --file_path="/hpc/dla_lti/dvanginneken/PLM-likelihoods/data/Kim/WU368_kim_bcr_heavy_elisa.csv"
echo "finished"

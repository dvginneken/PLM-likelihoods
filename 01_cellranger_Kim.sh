#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cellranger_Kim
#SBATCH --error=log/cellranger_Kim.error
#SBATCH --output=log/cellranger_Kim.out

export PATH=$PATH:/hpc/dla_lti/dvanginneken/cellranger-8.0.0/binå
cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/data/Kim/VDJ

cat SRR_Acc_List.txt | while read sample
do
    echo "start ${sample}"
    cellranger vdj --id=$sample \
         --reference=/hpc/dla_lti/dvanginneken/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0 \
         --fastqs=/hpc/dla_lti/dvanginneken/PLM-likelihoods/data/Kim/VDJ/fastq/${sample} \
         --sample=$sample
done

echo "finished"


#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=cellranger_Bruhn
#SBATCH --error=log/cellranger_Bruhn.error
#SBATCH --output=log/cellranger_Bruhn.out

cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/data/Bruhn/VDJ

cellranger vdj --id=Bruhn\
         --reference=/hpc/dla_lti/dvanginneken/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0 \
         --fastqs=/hpc/dla_lti/dvanginneken/PLM-likelihoods/data/Bruhn/fastq/SRR24140776 \
         --sample=SRR24140776

echo "finished"


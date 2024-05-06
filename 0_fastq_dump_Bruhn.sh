#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=FastqDump_Bruhn
#SBATCH --error=log/FastqDump_Bruhn.error
#SBATCH --output=log/FastqDump_Bruhn.out

cd /hpc/dla_lti/dvanginneken/PLM-likelihoods
#prefetch SRR24140776 --output-directory data/Bruhn/fastq/
fastq-dump --origfmt -I --split-files data/Bruhn/fastq/SRR24140776/SRR24140776.sra --outdir data/Bruhn/fastq/SRR24140776/


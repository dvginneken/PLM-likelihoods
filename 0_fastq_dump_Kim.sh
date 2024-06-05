#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=FastqDump_Kim
#SBATCH --error=log/FastqDump_Kim.error
#SBATCH --output=log/FastqDump_Kim.out

export PATH=$PATH:/hpc/dla_lti/dvanginneken/sratoolkit.3.1.0-centos_linux64/bin

cd /hpc/dla_lti/dvanginneken/PLM-likelihoods/data/Kim/VDJ
prefetch SRR17729674 --output-directory fastq/
fastq-dump --origfmt -I --split-files fastq/SRR17729674/SRR17729674.sra --outdir fastq/SRR17729674

#cat SRR_Acc_List.txt | while read srr
#do
#    prefetch $srr --output-directory fastq/
#    fastq-dump --origfmt -I --split-files fastq/${srr}/${srr}.sra --outdir fastq/${srr}
#done

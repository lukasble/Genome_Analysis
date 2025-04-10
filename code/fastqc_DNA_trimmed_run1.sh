#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J fastqc_dna_all
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<your_email@student.uu.se>
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/01_preprocessing/fastqc_trim/fastqc_dna_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/01_preprocessing/fastqc_trim/fastqc_dna_%j.err

# Load modules
module load bioinfo-tools
module load FastQC

# Ensure output directory exists
mkdir -p /home/lubl5753/Genome_Analysis/analyses/01_preprocessing/fastqc_trim

# Run FastQC on all .fastq.gz files in the DNA directory
fastqc -o /home/lubl5753/Genome_Analysis/analyses/01_preprocessing/fastqc_trim \
       /home/lubl5753/Genome_Analysis/data/raw_data/DNA_trimmed_linked/*.fastq.gz


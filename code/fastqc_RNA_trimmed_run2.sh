#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J fastqc_rna_trimmed
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/01_preprocessing/fastqc_trim/RNA_after_trim/fastqc_rna_trimmed_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/01_preprocessing/fastqc_trim/RNA_after_trim/fastqc_rna_trimmed_%j.err

# Load FastQC
module load bioinfo-tools
module load FastQC

# Set input/output directories
INPUT_DIR=/home/lubl5753/Genome_Analysis/data/trimmed_data
OUTPUT_DIR=/home/lubl5753/Genome_Analysis/analyses/01_preprocessing/fastqc_trim/RNA_after_trim

mkdir -p ${OUTPUT_DIR}

# Run FastQC on paired trimmed RNA files
fastqc -o ${OUTPUT_DIR} -t 2 \
    ${INPUT_DIR}/SRR4342137_1.paired.trimmed.fastq.gz \
    ${INPUT_DIR}/SRR4342137_2.paired.trimmed.fastq.gz \
    ${INPUT_DIR}/SRR4342139_1.paired.trimmed.fastq.gz \
    ${INPUT_DIR}/SRR4342139_2.paired.trimmed.fastq.gz


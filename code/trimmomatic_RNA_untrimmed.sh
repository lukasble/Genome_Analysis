#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J trim_rna
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/01_preprocessing/trimming_software/trim_rna_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/01_preprocessing/trimming_software/trim_rna_%j.err

# Load Trimmomatic
module load bioinfo-tools
module load trimmomatic

# Input and output directories
INPUT_DIR=/proj/uppmax2025-3-3/Genome_Analysis/3_Thrash_2017/RNA_untrimmed
OUTPUT_DIR=/home/lubl5753/Genome_Analysis/data/trimmed_data

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

echo "THIS IS THE UPDATED SCRIPT VERSION WITH A FOR LOOP"

# Loop through all *.1.fastq.gz files and find matching *.2.fastq.gz
for R1 in ${INPUT_DIR}/*.1.fastq.gz; do
    base=$(basename "$R1" .1.fastq.gz)
    R2=${INPUT_DIR}/${base}.2.fastq.gz

    echo "Processing $base"

    trimmomatic PE -threads 2 \
        ${R1} ${R2} \
        ${OUTPUT_DIR}/${base}_1.paired.trimmed.fastq.gz ${OUTPUT_DIR}/${base}_1.unpaired.fastq.gz \
        ${OUTPUT_DIR}/${base}_2.paired.trimmed.fastq.gz ${OUTPUT_DIR}/${base}_2.unpaired.fastq.gz \
        ILLUMINACLIP:/sw/apps/bioinfo/trimmomatic/0.39/rackham/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
done


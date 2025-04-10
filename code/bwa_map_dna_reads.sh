#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 08:00:00
#SBATCH -J bwa_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/slurm_log/bwa_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/slurm_log/bwa_%j.err

# Load modules
module load bioinfo-tools
module load bwa
module load samtools

# File paths
CONTIGS=/home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/megahit_out/final.contigs.fa
READS_DIR=/proj/uppmax2025-3-3/Genome_Analysis/3_Thrash_2017/DNA_trimmed
OUTDIR=/home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/mapping_out

mkdir -p ${OUTDIR}

# Index the contigs (only needs to be done once)
bwa index ${CONTIGS}

# Map and sort both DNA samples
bwa mem -t 4 ${CONTIGS} \
    ${READS_DIR}/SRR4342129_1.paired.trimmed.fastq.gz ${READS_DIR}/SRR4342129_2.paired.trimmed.fastq.gz |
    samtools sort -@ 4 -o ${OUTDIR}/SRR4342129.bam -

bwa mem -t 4 ${CONTIGS} \
    ${READS_DIR}/SRR4342133_1.paired.trimmed.fastq.gz ${READS_DIR}/SRR4342133_2.paired.trimmed.fastq.gz |
    samtools sort -@ 4 -o ${OUTDIR}/SRR4342133.bam -

# Merge BAM files (optional if doing co-assembly)
samtools merge ${OUTDIR}/merged.bam ${OUTDIR}/SRR4342129.bam ${OUTDIR}/SRR4342133.bam

# Index final BAM
samtools index ${OUTDIR}/merged.bam


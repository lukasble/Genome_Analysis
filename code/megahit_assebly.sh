#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 06:00:00
#SBATCH -J megahit_dna
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/slurm_log/megahit_dna_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/slurm_log/megahit_dna_%j.err

# Load module
module load bioinfo-tools
module load megahit

# Set input/output paths
READS_DIR=/proj/uppmax2025-3-3/Genome_Analysis/3_Thrash_2017/DNA_trimmed
OUTDIR=/home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/megahit_out

# Ensure SLURM log folder exists
mkdir -p /home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/slurm_log

# Run MEGAHIT
megahit \
  -1 ${READS_DIR}/SRR4342129_1.paired.trimmed.fastq.gz,${READS_DIR}/SRR4342133_1.paired.trimmed.fastq.gz \
  -2 ${READS_DIR}/SRR4342129_2.paired.trimmed.fastq.gz,${READS_DIR}/SRR4342133_2.paired.trimmed.fastq.gz \
  --k-list 41,61,81,101 \
  --kmin-1pass \
  --min-contig-len 1000 \
  -t 4 \
  -o ${OUTDIR}



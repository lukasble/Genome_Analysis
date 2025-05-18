#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J depth_summary
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out/depth_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out/depth_%j.err

# Load MetaBAT
module load bioinfo-tools
module load MetaBat/2.12.1

# Input & output
BAM=/home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/mapping_out/merged.bam
OUTDIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out

mkdir -p ${OUTDIR}

# Generate depth file
jgi_summarize_bam_contig_depths --outputDepth ${OUTDIR}/depth.txt ${BAM}


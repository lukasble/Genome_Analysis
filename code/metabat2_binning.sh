#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 02:00:00
#SBATCH -J metabat_binning
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out/metabat_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out/metabat_%j.err

# Load MetaBAT
module load bioinfo-tools
module load MetaBat/2.12.1

# Input
CONTIGS=/home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/megahit_out/final.contigs.fa
DEPTH=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out/depth.txt
OUTDIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out

# Run MetaBAT2
metabat2 -i ${CONTIGS} -a ${DEPTH} -o ${OUTDIR}/bin -t 4 --minContig 1500


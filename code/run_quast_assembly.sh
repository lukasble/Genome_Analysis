#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J quast_eval
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/quast_out/quast_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/quast_out/quast_%j.err

# Load module
module load bioinfo-tools
module load quast

# Input and output
CONTIGS=/home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/megahit_out/final.contigs.fa
OUTDIR=/home/lubl5753/Genome_Analysis/analyses/02_genome_assembly/quast_out

mkdir -p ${OUTDIR}

# Run QUAST
quast.py ${CONTIGS} -o ${OUTDIR} -t 2


#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J prokka_bins
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/prokka_out/prokka_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/prokka_out/prokka_%j.err

# Load modules
module load bioinfo-tools
module load prokka

# Define paths
BIN_DIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out
OUTDIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/prokka_out

# Make output directory
mkdir -p $OUTDIR

# List of bins to annotate (high-quality ones from CheckM)
bins=(
"bin.62.fa"
"bin.36.fa"
"bin.11.fa"
"bin.20.fa"
"bin.9.fa"
"bin.5.fa"
"bin.49.fa"
"bin.22.fa"
"bin.31.fa"
"bin.2.fa"
"bin.54.fa"
"bin.13.fa"
"bin.12.fa"
"bin.10.fa"
)

# Annotate each selected bin
for bin in "${bins[@]}"; do
    base=$(basename $bin .fa)
    echo "Annotating $base ..."
    
    prokka --outdir ${OUTDIR}/${base} \
           --prefix ${base} \
           --cpus 4 \
           --force \
           --compliant \
           ${BIN_DIR}/${bin}
done

echo "✅ Prokka annotation finished for all selected bins!"


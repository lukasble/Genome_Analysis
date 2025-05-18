#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -J gtdbtk_bins
#SBATCH --mem=128G               # Increased to 128 GB to prevent OOM
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/gtdbtk_out/gtdbtk_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/gtdbtk_out/gtdbtk_%j.err

# Load modules
module load bioinfo-tools
module load GTDB-Tk/2.4.0

# Define paths
BIN_DIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out
OUTDIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/gtdbtk_out

# Create output directory if it does not exist
mkdir -p ${OUTDIR}

# Run GTDB-Tk with optimized memory and CPU management
gtdbtk classify_wf \
    --genome_dir ${BIN_DIR} \
    --out_dir ${OUTDIR} \
    --cpus 16 \
    --pplacer_cpus 4 \
    --extension fa \
    --skip_ani_screen \
    --debug

echo "GTDB-Tk classification done!"


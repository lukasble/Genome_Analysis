#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 06:00:00
#SBATCH -J checkm_bins
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /proj/uppmax2025-3-3/nobackup/lubl5753/checkm_out/checkm_%j.out
#SBATCH -e /proj/uppmax2025-3-3/nobackup/lubl5753/checkm_out/checkm_%j.err

# Load modules
module load bioinfo-tools
module load CheckM/1.0.12

# Define paths
BIN_DIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out
OUTDIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/checkm_out

mkdir -p ${OUTDIR}

checkm lineage_wf \
  -x fa \
  -t 6 \
  -f ${OUTDIR}/checkm_summary.txt \
  --reduced_tree \
  ${BIN_DIR} \
  ${OUTDIR}

echo "CheckM analysis complete!"

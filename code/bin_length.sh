#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J bin_length_summary
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/05_abundance_analysis/abundance_out/bin_length_summary_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/05_abundance_analysis/abundance_out/bin_length_summary_%j.err

# Load necessary modules
module load bioinfo-tools

# Define paths
METABAT_DIR="/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out"
OUTDIR="/home/lubl5753/Genome_Analysis/analyses/05_abundance_analysis/abundance_out"
SUMMARY_FILE="${OUTDIR}/bin_length_summary.txt"

# Create output directory if not exists
mkdir -p ${OUTDIR}

# Initialize summary file
echo "Bin Total_Length" > ${SUMMARY_FILE}

# Loop through all bin files in Metabat directory
for BIN in ${METABAT_DIR}/*.fa; do
    BIN_NAME=$(basename ${BIN} .fa)
    TOTAL_LENGTH=$(grep -v ">" ${BIN} | wc -c)
    echo "$BIN_NAME $TOTAL_LENGTH" >> ${SUMMARY_FILE}
    echo "Processed $BIN_NAME with length $TOTAL_LENGTH"
done

echo "✅ Bin length calculation completed."


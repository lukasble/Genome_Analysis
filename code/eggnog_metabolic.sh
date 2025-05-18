#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8                # Use 8 cores for faster processing
#SBATCH -t 12:00:00         # Set sufficient time for all bins
#SBATCH -J eggnog_mapping
#SBATCH --mem=32G           # Use 32 GB for safe processing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/06_metabolic_pathways/eggnog_out/eggnog_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/06_metabolic_pathways/eggnog_out/eggnog_%j.err

# Load EggNOG Mapper (specific version)
module load bioinfo-tools
module load eggNOG-mapper/2.1.9

# Verify that EggNOG Mapper loaded correctly
if ! command -v emapper.py &> /dev/null; then
    echo "Error: EggNOG Mapper failed to load. Exiting."
    exit 1
fi

# Define paths
PROKKA_DIR="/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/prokka_out"
OUTDIR="/home/lubl5753/Genome_Analysis/analyses/06_metabolic_pathways/eggnog_out"

# Create output directory
mkdir -p ${OUTDIR}

# List of bins to process
bins=("bin.62.fa" "bin.36.fa" "bin.11.fa" "bin.20.fa" "bin.9.fa" "bin.5.fa" "bin.49.fa" "bin.22.fa" "bin.31.fa" "bin.2.fa" "bin.54.fa" "bin.13.fa" "bin.12.fa" "bin.10.fa")

# Summary file for metabolic pathways
SUMMARY_FILE="${OUTDIR}/eggnog_pathway_summary.txt"
echo "Bin KEGG_Pathway EC_Number Description" > $SUMMARY_FILE

# Loop through all bins
for bin in "${bins[@]}"; do
    base=$(basename $bin .fa)
    echo "Processing $base with EggNOG Mapper..."

    # Define paths
    PROKKA_FAA="${PROKKA_DIR}/${base}/${base}.faa"
    EGGNOG_OUTPUT="${OUTDIR}/${base}_eggnog"

    # Check if the protein file exists
    if [[ ! -f ${PROKKA_FAA} ]]; then
        echo "Error: Protein file ${PROKKA_FAA} not found. Skipping ${base}."
        continue
    fi

    # Run EggNOG Mapper
    emapper.py -i ${PROKKA_FAA} \
               -o ${EGGNOG_OUTPUT} \
               --cpu 8 \
               --output_dir ${OUTDIR} \
               --override  # Overwrites existing output

    # Parse EggNOG Mapper results for KEGG Pathways and EC Numbers
    ANNOTATIONS="${EGGNOG_OUTPUT}.annotations"
    if [[ -f ${ANNOTATIONS} ]]; then
        grep -v '^#' ${ANNOTATIONS} | awk -v bin="$base" '{if($7 != "-" && $8 != "-") print bin, $7, $8, $6}' >> ${SUMMARY_FILE}
    else
        echo "Error: Annotation file ${ANNOTATIONS} not found for ${base}." >> ${OUTDIR}/eggnog_errors.log
    fi
done

echo "✅ EggNOG Mapper analysis completed for all bins!"


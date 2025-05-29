#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J htseq_count
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH --output=/home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/htseq_out/htseq_%j.out
#SBATCH --error=/home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/htseq_out/htseq_%j.err

# Load modules
module load bioinfo-tools
module load htseq

# Define paths
GFF_DIR="/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/prokka_out2/gff_clean"
BAM_DIR="/home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/rna_mapping"
OUTPUT_DIR="/home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/htseq_out"
mkdir -p "$OUTPUT_DIR"

# Run HTSeq-count on each BAM file
for BAM_FILE in "${BAM_DIR}"/*_rna.bam; do
    BASENAME=$(basename "$BAM_FILE" _rna.bam)
    GFF_FILE="${GFF_DIR}/${BASENAME}.cleaned.gff"

    echo "→ Processing ${BASENAME} …"
    echo "Looking for GFF: ${GFF_FILE}"

    if [[ ! -f "$GFF_FILE" ]]; then
        echo "⚠️  GFF file not found: $GFF_FILE. Skipping $BASENAME."
        continue
    fi

    htseq-count \
        -f bam \
        -r pos \
        -s no \
        -t CDS \
        -i ID \
        "$BAM_FILE" \
        "$GFF_FILE" \
      > "${OUTPUT_DIR}/${BASENAME}_counts.txt"

    echo "✔️  Counts written to ${OUTPUT_DIR}/${BASENAME}_counts.txt"
done

# -----------------------------------------------------
# Summary of total read counts per sample
# -----------------------------------------------------
echo ""
echo "===== Summary of Read Counts Per Sample ====="
cd "$OUTPUT_DIR"
for f in *_counts.txt; do
    TOTAL=$(grep -v '^__' "$f" | awk '{sum += $2} END {print sum}')
    echo "$(basename "$f"): $TOTAL reads counted (excluding special entries)"
done


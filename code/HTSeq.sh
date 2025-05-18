#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J htseq_counts
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/htseq_out/htseq_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/htseq_out/htseq_%j.err

module load bioinfo-tools
module load htseq/2.0.2

BAM_DIR=/home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/rna_mapping
PROKKA_DIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/prokka_out
OUT_DIR=/home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/htseq_out

mkdir -p ${OUT_DIR}

BINS=("bin.62" "bin.36" "bin.11" "bin.20" "bin.9" "bin.5" "bin.49" "bin.22" "bin.31" "bin.2" "bin.54" "bin.13" "bin.12" "bin.10")

for BIN in "${BINS[@]}"; do
    GFF_FILE=${PROKKA_DIR}/${BIN}/${BIN}.gff
    BAM_FILE=${BAM_DIR}/${BIN}_rna.bam
    OUT_FILE=${OUT_DIR}/${BIN}_counts.txt

    if [[ -f "$GFF_FILE" && -f "$BAM_FILE" ]]; then
        # Run HTSeq-count for CDS only
        htseq-count -f bam -r pos -i ID "$BAM_FILE" <(grep -v '^#' "$GFF_FILE" | awk '$3 == "CDS"') > "$OUT_FILE"
    else
        echo "Skipping $BIN: missing GFF or BAM"
    fi
done

echo "HTSeq-count processing for CDS only completed."


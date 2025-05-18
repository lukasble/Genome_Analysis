#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J rna_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/rna_mapping/rna_mapping_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/rna_mapping/rna_mapping_%j.err

# Load required modules
module load bioinfo-tools
module load bwa
module load samtools

# Define input and output paths
BIN_DIR="/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out"
READ1="/home/lubl5753/Genome_Analysis/data/trimmed_data/SRR4342137_1.paired.trimmed.fastq.gz"
READ2="/home/lubl5753/Genome_Analysis/data/trimmed_data/SRR4342137_2.paired.trimmed.fastq.gz"
OUTDIR="/home/lubl5753/Genome_Analysis/analyses/04_functional_annotation/rna_mapping"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Loop through each bin and map RNA reads
for BIN in ${BIN_DIR}/bin.*.fa; do
    BIN_NAME=$(basename "$BIN" .fa)
    echo "Processing $BIN_NAME..."

    # Index the bin
    bwa index "$BIN"

    # Align reads and sort output BAM
    bwa mem -t 4 "$BIN" "$READ1" "$READ2" | \
        samtools view -bS - | \
        samtools sort -@ 4 -o "${OUTDIR}/${BIN_NAME}_rna.bam" -

    # Index the BAM file
    samtools index "${OUTDIR}/${BIN_NAME}_rna.bam"
done

echo "RNA mapping complete!"


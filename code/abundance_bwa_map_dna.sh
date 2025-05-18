#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8                # Use 8 cores (4 for mapping, 4 for sorting)
#SBATCH -t 12:00:00         # Increased time for multiple bins
#SBATCH -J abundance_analysis_v2
#SBATCH --mem=32G           # Increased memory for safe sorting
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /crex/proj/uppmax2025-2-288/nobackup/lubl5753/nextflow_lab/abudance_analysis/abundance_%j.out
#SBATCH -e /crex/proj/uppmax2025-2-288/nobackup/lubl5753/nextflow_lab/abudance_analysis/abundance_%j.err

# Load required modules
module load bioinfo-tools
module load bwa
module load samtools

# Define paths
BIN_DIR="/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out"
READS_DIR="/home/lubl5753/Genome_Analysis/data/raw_data/DNA_trimmed_linked"
OUTDIR="/crex/proj/uppmax2025-2-288/nobackup/lubl5753/nextflow_lab/abudance_analysis"

# Create output directory
mkdir -p ${OUTDIR}

# List of bins to process
bins=("bin.62.fa" "bin.36.fa" "bin.11.fa" "bin.20.fa" "bin.9.fa" "bin.5.fa" "bin.49.fa" "bin.22.fa" "bin.31.fa" "bin.2.fa" "bin.54.fa" "bin.13.fa" "bin.12.fa" "bin.10.fa")

# List of samples to process
samples=("SRR4342129" "SRR4342133")

# Summary file
SUMMARY_FILE="${OUTDIR}/abundance_summary.txt"
echo "Sample Bin Total_Reads Mapped_Reads Mapped_Percentage" > $SUMMARY_FILE

# Loop over samples and bins
for sample in "${samples[@]}"; do
    for bin in "${bins[@]}"; do
        base=$(basename $bin .fa)
        echo "Processing $sample with $base ..."

        # Define output BAM file
        BAM_FILE="${OUTDIR}/${sample}_${base}.sorted.bam"

        # Mapping with BWA
        bwa mem -t 4 ${BIN_DIR}/${bin} ${READS_DIR}/${sample}_1.paired.trimmed.fastq.gz ${READS_DIR}/${sample}_2.paired.trimmed.fastq.gz | \
        samtools sort -m 4G -@ 4 -o ${BAM_FILE} -T ${OUTDIR}/${sample}_${base}_temp

        # Calculating mapping statistics
        total_reads=$(samtools view -c ${BAM_FILE})
        mapped_reads=$(samtools view -c -F 4 ${BAM_FILE})
        mapped_percent=$(echo "scale=2; (${mapped_reads}/${total_reads})*100" | bc)

        # Adding to summary file
        echo "${sample} ${base} ${total_reads} ${mapped_reads} ${mapped_percent}" >> ${SUMMARY_FILE}

        # Cleanup temp files
        rm -f ${OUTDIR}/${sample}_${base}_temp*
    done
done

echo "✅ Abundance analysis completed for all bins and samples!"

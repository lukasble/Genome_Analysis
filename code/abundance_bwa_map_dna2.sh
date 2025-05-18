#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J abundance_analysis
#SBATCH --mem=32G  # Increased memory to 32G for safe sorting
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/05_abundance_analysis/abundance_out/abundance_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/05_abundance_analysis/abundance_out/abundance_%j.err

# Load required modules
module load bioinfo-tools
module load bwa
module load samtools

# Define paths
BIN_PATH="/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out/bin.62.fa"
READS_DIR="/home/lubl5753/Genome_Analysis/data/raw_data/DNA_trimmed_linked"
OUTDIR="/home/lubl5753/Genome_Analysis/analyses/05_abundance_analysis/abundance_out"

mkdir -p ${OUTDIR}

# Clean any existing corrupted BAM files
rm -f ${OUTDIR}/test.sorted.bam ${OUTDIR}/test.resorted.bam ${OUTDIR}/test.bam

echo "Mapping reads to bin.62 with BWA..."
bwa mem -t 4 $BIN_PATH ${READS_DIR}/SRR4342129_1.paired.trimmed.fastq.gz ${READS_DIR}/SRR4342129_2.paired.trimmed.fastq.gz 2> ${OUTDIR}/test_bwa.log | \
samtools view -Sb - > ${OUTDIR}/test.bam

echo "Sorting BAM file..."
samtools sort -m 4G -@ 4 -o ${OUTDIR}/test.sorted.bam ${OUTDIR}/test.bam 2> ${OUTDIR}/test_sort.log
rm ${OUTDIR}/test.bam  # Clean up to save space

echo "Verifying BAM file integrity..."
samtools quickcheck -v ${OUTDIR}/test.sorted.bam > ${OUTDIR}/test_check.log 2>&1

# If the BAM file is corrupted, delete and re-map
if [[ -s ${OUTDIR}/test_check.log ]]; then
    echo "BAM file is corrupted. Restarting mapping and sorting from scratch..."
    rm -f ${OUTDIR}/test.sorted.bam
    bwa mem -t 4 $BIN_PATH ${READS_DIR}/SRR4342129_1.paired.trimmed.fastq.gz ${READS_DIR}/SRR4342129_2.paired.trimmed.fastq.gz 2> ${OUTDIR}/test_bwa.log | \
    samtools view -Sb - > ${OUTDIR}/test.bam
    samtools sort -m 4G -@ 4 -o ${OUTDIR}/test.sorted.bam ${OUTDIR}/test.bam
    rm ${OUTDIR}/test.bam
else
    echo "BAM file is intact."
fi

echo "Generating flagstat..."
samtools flagstat ${OUTDIR}/test.sorted.bam > ${OUTDIR}/test.flagstat.txt

echo "Calculating depth..."
samtools depth ${OUTDIR}/test.sorted.bam > ${OUTDIR}/test.depth.txt

echo "Final Results:"
cat ${OUTDIR}/test.flagstat.txt
head ${OUTDIR}/test.depth.txt


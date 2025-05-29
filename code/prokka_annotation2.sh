#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J prokka_bins2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lukas.bleichner.5753@student.uu.se
#SBATCH -o /home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/prokka_out2/prokka_%j.out
#SBATCH -e /home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/prokka_out2/prokka_%j.err

module load bioinfo-tools
module load prokka

# input FASTAs
BIN_DIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/metabat_out

# fresh output folder
OUTDIR=/home/lubl5753/Genome_Analysis/analyses/03_structural_annotation/prokka_out2
mkdir -p "$OUTDIR"

# only these 14 bins:
bins=(
  bin.62 bin.36 bin.11 bin.20 bin.9 bin.5
  bin.49 bin.22 bin.31 bin.2 bin.54 bin.13
  bin.12 bin.10
)

for b in "${bins[@]}"; do
  echo "🚀 Annotating $b..."
  prokka \
    --outdir "${OUTDIR}/${b}" \
    --prefix "${b}" \
    --cpus 4 \
    "${BIN_DIR}/${b}.fa"

  # now clean the GFF (strip off FASTA, ready for HTSeq)
  awk '/^##FASTA/{exit} {print}' \
      "${OUTDIR}/${b}/${b}.gff" \
    > "${OUTDIR}/${b}/${b}.cleaned.gff"

  echo "✅ $b: cleaned GFF written to ${b}.cleaned.gff"
done

echo "🎉 All done.  Prokka + GFF‐clean for HTSeq ready in prokka_out2/"


# Genome Analysis: Northern Gulf of Mexico

This repository contains Bash/SLURM scripts for a metagenomic and metatranscriptomic analysis inspired by the paper [“Metabolic roles of uncultivated bacterioplankton lineages in the Northern Gulf of Mexico ‘dead zone’”](https://pmc.ncbi.nlm.nih.gov/articles/PMC5596340/).

The scripts were written to run on UPPMAX HPC at Uppsala University. They document an attempted workflow from read quality control through assembly, binning, genome annotation, RNA read counting, and abundance analysis. The repository contains scripts and a progress log, but no input data, result files, or final conclusions; completion of the analysis therefore cannot be confirmed from this repository alone.

## Project Contents

- [`code/`](code/) contains the analysis scripts. Most are SLURM batch jobs and load software with the UPPMAX `module` system.
- [`README_log.md`](README_log.md) is the dated work log.
- `data/` and `analyses/` are excluded from version control (see `.gitignore`). The input data and generated outputs are not included here.

## Data Referenced by the Scripts

The scripts refer to paired-end SRA read sets:

| Data | Accession(s) | Use shown in scripts |
| --- | --- | --- |
| DNA | `SRR4342129`, `SRR4342133` | Co-assembly, DNA mapping, and bin abundance estimates |
| RNA | `SRR4342137`, `SRR4342139` | RNA trimming and trimmed-read quality control |
| RNA | `SRR4342137` | Mapping to selected bins and HTSeq counting |

These roles are based on script inputs, not an independently verified sample sheet. Confirm the sample metadata and experimental design before interpreting or reusing the analysis.

## Workflow Represented

1. **Read quality control:** FastQC on trimmed DNA and untrimmed RNA reads.
2. **RNA trimming:** Trimmomatic paired-end trimming with adapter clipping, quality trimming, and a 50 bp minimum length; FastQC is then run on trimmed reads.
3. **Assembly:** MEGAHIT co-assembly of the two DNA read sets, retaining contigs of at least 1,000 bp. QUAST evaluates the assembly.
4. **Mapping and binning:** BWA-MEM maps DNA reads to the assembly, samtools prepares and merges BAM files, and MetaBAT2 bins contigs using summarized depth information. The MetaBAT2 command uses a 1,500 bp minimum contig length.
5. **Bin quality and annotation:** CheckM evaluates bin completeness/contamination; Prokka annotates a selected set of 14 bins; GTDB-Tk classifies bins taxonomically.
6. **RNA analysis:** BWA maps one trimmed RNA read set to the selected bins, and HTSeq-count is configured to count CDS features from Prokka GFF files.
7. **Functional and abundance analyses:** An EggNOG-mapper script attempts to summarize KEGG pathway and EC annotations. A separate BWA/samtools script estimates DNA read mapping percentages per bin and sample.


## Project Wiki

See the separately maintained project GitHub wiki for the detailed guide to connecting to UPPMAX from Linux and running the analysis. This README is a concise overview of the repository, workflow, and scripts. The wiki can be found above.

## Script Index

| Script | Purpose |
| --- | --- |
| `fastqc_dna_trimmed.sh` | FastQC for trimmed DNA reads |
| `fastqc_rna_untrimmed.sh` | FastQC for untrimmed RNA reads |
| `trimmomatic_RNA_untrimmed.sh` | Paired-end trimming of RNA reads |
| `fastqc_rna_trimmed.sh` | FastQC after RNA trimming |
| `megahit_assebly.sh` | MEGAHIT DNA co-assembly |
| `run_quast_assembly.sh` | QUAST assembly assessment |
| `bwa_map_dna_reads.sh` | DNA read mapping and BAM merge |
| `generate_depth_file.sh` | Contig depth summary for binning |
| `metabat2_binning.sh` | MetaBAT2 binning |
| `checkm_quality.sh` | CheckM bin quality assessment |
| `bin_length.sh` | Calculate bin sequence lengths |
| `prokka_annotation2.sh` | Prokka annotation of selected bins and cleaned GFF generation |
| `GTDB-tk.sh` | GTDB-Tk taxonomic classification |
| `bwa_map_rna_reads.sh` | RNA read mapping to selected bins |
| `HTSeq.sh` | HTSeq feature counts from mapped RNA reads |
| `eggnog_metabolic.sh` | EggNOG-mapper annotations and pathway summary attempt |
| `abundance_bwa_map_dna.sh` | DNA mapping percentages per sample and bin |

## Known Reproducibility Issues

- Paths vary among scripts and refer to multiple project/storage locations. They need to be reconciled with the actual UPPMAX data layout.
- The Prokka script writes to `prokka_out2` and creates cleaned GFFs inside each bin's output directory. The HTSeq script expects cleaned GFFs in a separate `prokka_out2/gff_clean` directory. Confirm or correct this handoff before running HTSeq.
- The EggNOG script reads proteins from `prokka_out`, whereas the available Prokka script writes to `prokka_out2`. Confirm the intended annotation directory and that the expected protein files exist.
- The abundance script uses a different UPPMAX project path from most scripts, and its reads path is named `DNA_trimmed_linked` while the filenames requested are trimmed-paired FASTQs. Verify the data there before running it.
- The scripts do not include a unified configuration file, environment lock, sample sheet, or captured software/database versions for every tool. Record these, along with job outcomes and key result files, for a reproducible final analysis.

## Current Status

The analysis results are documented in a separate report and are not included in this repository. They may be published on GitHub in the future.
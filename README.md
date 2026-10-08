# Genome Analysis: Northern Gulf of Mexico

![Linux](https://img.shields.io/badge/OS-Linux-blue?logo=linux)
![HPC](https://img.shields.io/badge/HPC-UPPMAX%20%28Rackham%29-red)
![Workload Manager](https://img.shields.io/badge/Scheduler-SLURM-orange)
![Domain](https://img.shields.io/badge/Domain-Metagenomics%20%26%20Metatranscriptomics-green)

> **Key Tech Stack:** `Linux` • `SLURM` • `UPPMAX HPC (Uppsala University)` • `Bash` • `Bioinformatics Pipelines`

This repository contains Bash/SLURM scripts for a metagenomic and metatranscriptomic analysis executed on the **UPPMAX HPC (Rackham cluster)** at Uppsala University. The workflow is inspired by the research paper [“Metabolic roles of uncultivated bacterioplankton lineages in the Northern Gulf of Mexico ‘dead zone’”](https://pmc.ncbi.nlm.nih.gov/articles/PMC5596340/).

---

## Compute Environment & Tools

* **Infrastructure & OS:** Remote HPC via Linux / SSH connection to UPPMAX (Uppsala Multidisciplinary Center for Advanced Computational Science)
* **Workload Management:** SLURM batch scripts (`#SBATCH` directives, parallel execution, resource requests)
* **Software Management:** UPPMAX `module` system
* **Bioinformatics Suite:**
  * **QC & Trimming:** FastQC, Trimmomatic
  * **Assembly & Assessment:** MEGAHIT, QUAST
  * **Mapping & Binning:** BWA-MEM, Samtools, MetaBAT2
  * **Quality & Annotation:** CheckM, Prokka, GTDB-Tk, EggNOG-mapper
  * **Expression:** HTSeq-count

---

## Project Contents

- [`code/`](code/) contains the analysis scripts. Most are SLURM batch jobs loading software via UPPMAX modules.
- [`README_log.md`](README_log.md) is the dated work log detailing the HPC execution timeline.
- `data/` and `analyses/` are excluded from version control (see `.gitignore`). Input datasets and generated outputs are omitted due to storage constraints.

## Data Referenced by the Scripts

The scripts process paired-end SRA read sets:

| Data | Accession(s) | Use shown in scripts |
| --- | --- | --- |
| DNA | `SRR4342129`, `SRR4342133` | Co-assembly, DNA mapping, and bin abundance estimates |
| RNA | `SRR4342137`, `SRR4342139` | RNA trimming and trimmed-read quality control |
| RNA | `SRR4342137` | Mapping to selected bins and HTSeq counting |

> [!NOTE]
> These scripts were developed for academic research workflows on UPPMAX. Hardcoded absolute pathing or module versions may require adjustment for execution on other HPC environments.

---

## Workflow Represented

1. **Read Quality Control:** FastQC on trimmed DNA and untrimmed RNA reads.
2. **RNA Trimming:** Trimmomatic paired-end trimming with adapter clipping and 50 bp minimum length filter; FastQC post-trimming.
3. **Co-Assembly:** MEGAHIT co-assembly of DNA read sets (contigs ≥ 1,000 bp), evaluated via QUAST.
4. **Mapping & Binning:** BWA-MEM mapping to assembly, Samtools BAM processing, and MetaBAT2 contig binning using coverage profiles.
5. **Bin Assessment & Annotation:** CheckM completeness/contamination assessment; Prokka annotation on selected bins; GTDB-Tk taxonomic assignment.
6. **Metatranscriptomics (RNA):** BWA mapping of trimmed RNA reads to selected bins and HTSeq feature counting against Prokka GFFs.
7. **Functional & Abundance Analysis:** EggNOG-mapper pathway annotation (KEGG/EC) and BWA/Samtools relative abundance estimation.

---

## Project Wiki

See the separately maintained **GitHub Wiki** above for a detailed guide on establishing SSH connections to UPPMAX from Linux and orchestrating the SLURM job execution.

---

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

---

## Current Status

The analysis results are documented in a separate report and are omitted from this repository due to size constraints.

The analysis results are documented in a separate report and are omitted from this repository due to size constraints.

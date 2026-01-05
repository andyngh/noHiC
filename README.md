# noHiC
Pangenome-based contig scaffolding pipeline with no HiC data 

## 1. Introduction

### 1.1. noHiC workflow

noHiC is a reference-guided genome assembly pipeline covering essential steps in contig scaffolding with four independent sub-scripts, including contig decontamination (`nohic-clean`), misassembly correction, reference-guided contig scaffolding (`nohic-refpick` and `nohic-asm`), and assembly evaluation (`nohic-eval`). Given an input FASTA file containing contig-level assembly, the pipeline starts with removing contigs not originating from a specified taxonomic group and optionally organellar contigs. After that, misassemblies in the selected contigs will be corrected in three different approaches using **error-corrected** longs reads (each can be skipped), including clip-based chimeric contig breakings by [CRAQ](https://github.com/JiaoLaboratory/CRAQ/tree/main), small misassembly (base substitutions, expansions, collapses, and haplotype switches) corrections by [Inspector](https://github.com/Maggi-Chen/Inspector), and reference-guided misjoin breakings by [RagTag correct](https://github.com/malonge/RagTag/wiki/correct). The corrected contigs are then scaffolded based on a reference genome with [RagTag scaffold](https://github.com/malonge/RagTag/wiki/scaffold). The most important feature of noHiC is that users can utilized publicly available pangenome graphs built by tools like [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) to synthesized a genetically closed, personalized reference genome (synref) for reference-guided contig error correction and scaffolding based on [KMC](https://github.com/refresh-bio/KMC?tab=readme-ov-file) and [vg haplotype](https://github.com/vgteam/vg/wiki/Haplotype-Sampling), minimizing the amount of false contig correction. Quality of the scaffolded assembly can be assessed based on both metrics (N50, auN, gap number, BUSCO...) and visualizations (i.e. comparing target assembly with a reference genome using a dot plot built by [paf2dotplot](https://github.com/moold/paf2dotplot) and visualizing the positions of remaining misassemblies on chromosomes by `nohic-viz.R`).   

![noHiC workflow](https://github.com/andyngh/noHiC/blob/main/Fig1.png)

### 1.2. dependency citations

### 1.3. noHiC citation

## 2. Installation

## 3. Tutorial

### 3.1. nohic-clean.sh: Removing Contaminant Contigs

General usage:

```
nohic-clean.sh -i <contig.fa> -a <adapter.txt> -o <common_out> -t <threads> -kp <kp_prefix> -tg <taxonomic_group> [-ros <mt.cl.fasta>] [-ioc yes|no] [-m yes|no] [-kdb <kraken2_db_dir>] [--resume]
Required:
  -i, --input-fasta <.fa|.fasta>                                FASTA file with raw input contigs
  -a, --adapter-sequence <.txt>                                 A text file containing adapter sequences (one adapter per line)
  -o, --output-directory <dir>                                  Common output directory for storing results from contig cleaning steps
  -t, --threads <int>                                           Thread number for kraken2 and blastn
  -kp, --prefix-for-taxonomic-classification <str>              Prefix for kraken2 outputs (e.g. sample_1)
  -tg, --taxonomic-group <str>                                  Include contigs from this taxonomic group (e.g. Viridiplantae, Bacteria, Mammalia...)
Optional:
  -ros, --reference-organellar-sequences <.fa|.fasta>           FASTA file containing reference organellar DNA sequences
  -ioc, --identify-organellar-contigs <yes|no>                  Check if there is organellar DNA in the contigs (Default: no)
  -m,  --memory-mapping <yes|no>                                Use --memory-mapping option of kraken2 (Default: yes)
  -kdb, --kraken2-database <dir>                                Path to the directory containing a kraken2 database (*.k2d files) (required if -m no)
  --resume                                                      Resume the pipeline from the earliest failed steps
  -h|--help                                                     Display this help message
```
**Purposes**

The initial steps of the noHic pipeline involve removal of contigs originating from non-target species (contaminant contigs) and organelles, thereby making sure that these contigs will not intefere with the downstream scaffolding. The nohic-clean first checks for the presence of adapters in the input contigs. If there are adapters, the script will stop with an error message. It is highly recommended that users should trim their reads with some effective programs like [TGSFilter](https://github.com/HuiyangYu/TGSFilter) and re-assemble their contigs before executing nohic-clean again with the `--resume` flag. Once identifying that the contigs are adapter-free, nohic-clean calls Kraken2 and TaxonKit to determine contaminant contigs in the input. If required by users, nohic-clean will also identify contigs from organelles based on BLASTn and a FASTA file containing reference mtDNA and cpDNA sequences. The contigs marked as contaminants will finally be removed.       

**Adapter detection**



**Contig taxonomic classification options**



**Organellar DNA removal options**



**Outputs** 



### 3.2. nohic-refpick.sh: Creating Personalized Reference for the Cleaned Contigs

### 3.3. nohic-asm.sh: Scaffolding Cleaned Contigs Based on a Reference Genome

### 3.4. nohic-eval.sh: Assembly Evaluation and Visualization

### 3.5. Example - Scaffolding the *A. thaliana* CAMA-C-2 Genome

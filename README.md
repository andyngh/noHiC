# noHiC
Pangenome-based contig scaffolding pipeline with no HiC data 

## 1. Introduction

### 1.1. noHiC workflow

noHiC is a reference-guided genome assembly pipeline covering essential steps in contig scaffolding with four independent sub-scripts, including contig decontamination (`nohic-clean`), misassembly correction, reference-guided contig scaffolding (`nohic-refpick` and `nohic-asm`), and assembly evaluation (`nohic-eval`). Given an input FASTA file containing contig-level assembly, the pipeline starts with removing contigs not originating from a specified taxonomic group and optionally organellar contigs. After that, misassemblies in the selected contigs will be corrected in three different approaches using **error-corrected** longs reads (each can be skipped), including clip-based chimeric contig breakings by [CRAQ](https://github.com/JiaoLaboratory/CRAQ/tree/main), small misassembly (base substitutions, expansions, collapses, and haplotype switches) corrections by [Inspector](https://github.com/Maggi-Chen/Inspector), and reference-guided misjoin breakings by [RagTag correct](https://github.com/malonge/RagTag/wiki/correct). The corrected contigs are then scaffolded based on a reference genome with [RagTag scaffold](https://github.com/malonge/RagTag/wiki/scaffold). The most important feature of noHiC is that users can utilized publicly available pangenome graphs built by tools like [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) to synthesized a genetically closed, personalized reference genome (synref) for reference-guided contig error correction and scaffolding based on [KMC](https://github.com/refresh-bio/KMC?tab=readme-ov-file) and [vg haplotype](https://github.com/vgteam/vg/wiki/Haplotype-Sampling), minimizing the amount of false contig correction. Quality of the scaffolded assembly can be assessed based on both metrics (N50, auN, gap number, BUSCO...) and visualizations (i.e. comparing target assembly with a reference genome using a dot plot built by [paf2dotplot](https://github.com/moold/paf2dotplot) and visualizing the positions of remaining misassemblies on chromosomes by `nohic-viz.R`). Please find the detailed descriptions of noHiC in **our paper**.  

![noHiC workflow](https://github.com/andyngh/noHiC/blob/main/Fig1.png)

### 1.2. dependency citations

### 1.3. noHiC citation

## 2. Installation

## 3. Tutorial

To guide users through the sub-scripts of noHiC, we have provided a small example of *A. thaliana* CAMA-C-2 contig scaffolding, with the essential files in https://github.com/andyngh/noHiC/blob/main/example. In this example, we will try to replicate the manually curated CAMA-C-2 assembly from [this paper](https://doi.org/10.1038/s41586-023-06062-z) using noHiC.

### 3.1. Preparations of HiFi Reads

The FASTQ file containing HiFi reads of *A. thaliana* CAMA-C-2 is downloaded from NCBI SRA using `prefetch` and `fastq-dump` from [SRA Toolkit](https://github.com/ncbi/sra-tools). We highly recommend [TGSFilter](https://github.com/HuiyangYu/TGSFilter) for trimming adapters in the reads. 

```
prefetch --max-size 200G ERR10084604

fastq-dump --origfmt ./ERR10084604

tgsfilter -i ERR10084604.fastq -o CAMA-C-2-hifi_reads.ALL.trimmed.fastq.gz -x hifi -t 24
```

In the benchmarks of **our paper**, 90% of the CAMA-C-2 reads will be used for assembly and 10% will be used **only** for assembly evaluation. Though **you do not need to do this step in your own assembly project**, we will conduct the read subsampling here so that the tutorial's results will be the same as in our paper.   

```
# 10% for evaluation
$nohic/seqkit sample -p 0.10 -j 100 -s 3108 CAMA-C-2-hifi_reads.ALL.trimmed.fastq.gz -o CAMA-C-2-hifi_reads.EVAL.trimmed.fastq.gz

# 90% for assembly
$nohic/seqkit grep -v -f <(zcat CAMA-C-2-hifi_reads.EVAL.trimmed.fastq.gz | sed -n '1~4s/^@//p') \
  CAMA-C-2-hifi_reads.ALL.trimmed.fastq.gz -o CAMA-C-2-hifi_reads.FOR_ASM.trimmed.fastq.gz
```

### 3.2. Contaminant Contig Removal

The `nohic-clean` sub-script will be used to remove contaminant contigs from an input contig-level assembly. The general usage of the script is as following:

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

Before running `nohic-clean`, we need to download the `core-nt` database for contig taxonomic classifications. 

```
wget https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz

tar -xzf ./k2_core_nt_20251015.tar.gz
```

Now we remove any contigs that are not from green plants (Viridiplantae) and organellar contigs.

```
nohic-clean.sh -i CAMA-C-2.asm.bp.p_ctg.fa -a PacBio_adapter.txt -o CAMA-C-2.asm.bp.p_ctg.cleaning -t 30 -kp CAMA-C-2 -tg Viridiplantae -ioc yes -ros mt.cl.fasta \
               -m no -kdb absolute/path/to/kraken2_core_nt_db/
```

Where, `CAMA-C-2.asm.bp.p_ctg.fa` is the input contigs. `PacBio_adapter.txt` contains adapter sequences with no headers and one sequence per line. The reference mt- and cpDNA of *A. thaliana* are in `mt.cl.fasta`.

**Note:** If you have root priviledge on your machine and many contig-level assemblies to decontaminate, we recommend executing the following commands before running `nohic-clean` with `-m yes` and **no** `-kdb`.

```
# Resize your /dev/shm directory to fit the core_nt database

# For example, if your kraken2 database size is 242G, run this:

sudo mount -o remount,size=244G /dev/shm

# Copy your *.k2d files to /dev/shm/

cp /path/to/your/database/directory/*.k2d /dev/shm/
```

**Outputs**

`nohic-clean.sh` will create four sub-directories inside the common output directory `CAMA-C-2.asm.bp.p_ctg.cleaning`. The output files in the subdirectories are as following.

*1_adapter_content_check*

In the case where there is no adapter in the contigs:

`adapter_check.log`: the log file of the adapter detection step.

`adapter_content.txt`: this file shows adapter detection result. It will say "None" if there is no adapter.

`fgrep_matches.txt`: this file is empty in this case.

`step_1_done.txt`: step completion marker. It will say "ok" if there is no adapter detected.

In the case where there are adapters in the contigs:

The sub-script will stop if it detects adapters.

`adapter_check.log`: the log file of the adapter detection step.

`adapter_matches.fa` and `adapter_matches.txt`: a FASTA file containing the contigs with adapters and a text file with the names of adapter-containing contigs, respectively. Users can choose to remove these contigs out of the input assembly before re-running `nohic-clean` if they are small and/or there are only a few of them. If there are long and/or many adapter-containing contigs, we highly recommend you to do read adapter trimming again with a more effective tool and re-assembly the contigs before running `nohic-clean`.

`fgrep_matches.txt`: this file now tell you where the adapters are detected in each contig in the format: input_file_path:in-contig_line_number_of_detected_adapter:contig_sequence

`step_1_done.txt`: this file will say "failed" in the case of adapters found.

*2_contig_taxonomic_classification*

`CAMA-C-2.kr` and `CAMA-C-2.report`: output files from Kraken2-based contig taxonomic classifications. Please see this [link](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats) to understand the output formats.

`CAMA-C-2.lineage`: this file is from Taxonkit. It has lineage information created using the TaxIDs of contigs, identified by Kraken2. 

`name_of_contaminant_contigs.txt`: this file contains the names of contaminant contigs. It will be used for decontamination later.

`kraken2_and_taxonkit.log`: the log file of the step.

`step_2_done.txt`: step completion marker

*3_organellar_contigs_identification*

`blast_out.tab`: this file contains unfiltered BLAST results with contigs as queries and mt- and cpDNA sequences as subjects in [output format 6](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6).

`filter_blast_out.tab`: this file contains the filtered BLAST results with pident and qcovs of >= 90%. 

`filter_blast.py`: the helper script used to filter the BLAST results.

`organellar_dna_containing_contigs.txt`: this file contains the names of organellar contigs. It will be used for decontamination later.

`blast_formatdb_and_filter.log`: log file of the step

`step_3_done.txt`: step completion marker

*4_contig_purification*

`CAMA-C-2.asm.bp.p_ctg.pure.fa`: the FASTA files with contaminant contigs removed. This is the **main output file**.

`contigs_from_contaminant_n_organelle.txt`: A list combining `name_of_contaminant_contigs.txt` in step 2 and `organellar_dna_containing_contigs.txt` in step 3.

`remove_contamination.log`: the log file of the step

`remove_contamination.py`: the helper script used to remove contaminant contigs.

`step_4_done.txt`: step completion marker

**Note:** If any of the step of `nohic-clean` failed, you can resume the sub-script with the same command plus the `--resume` flag after fixing the errors and deleting every file in the sub-directory of the previous failed step. 

To prepare the clean contigs for scaffolding, we will exclude contigs shorter than 100 kb (as mentioned in [this paper](https://doi.org/10.1038/s41586-023-06062-z)).

```
seqkit seq -m 100000 CAMA-C-2.asm.bp.p_ctg.pure.fa -o CAMA-C-2.asm.bp.p_ctg.pure.for_asm.fa
```

### 3.2. Creating Personalized Reference for the Target Genome

In parallel with contig decontamination, we will create a personalized reference (synref) for *A. thaliana* CAMA-C-2 here using `nohic-refpick.sh`. The general usage of the subscript is as following.

```
nohic-refpick.sh -s <hifi_reads.fastq.gz|target_contigs.fa> -g <pangenome_graph.gbz> -i <pangenome_graph.hapl> -o <outprefix> [-t <threads>] [-m <memory>] [-v <0|1|2|3>] [-p <yes|no>] [-r <reference.fa>]
Required:
  -s, --input-sequence <.fastq/.fasta>  Input contig assembly or fastq file with long reads (fastq recommended)
  -g, --gbz     <.gbz>                  Input pangenome graph in GBZ format
  -i, --hapidx  <.hapl>                 Haplotype information from the pangenome graph
  -o, --outprefix <str>                 Prefix for outputs

Optional:
  -t, --threads <int>                   Thread number (default: 1)
  -m, --ram-gb  <int>                   Maximum memory for KMC in GB (default: 32)
  -k, --kmer    <int>                   k-mer size for KMC (default: 29)
  -v, --verbosity <int>                 vg verbosity level (0 = silent, 1 = basic, 2 = detailed, 3 = debug; default: 2)
  -p, --patch <yes|no>                  Run personalized reference patching step (default: yes)
  -r, --reference <.fasta>              Reference fasta for patching (required if --patch yes)
  -h, --help                            Display this help message
      --version                         Display version number
```

The following command will be run to create a synref for CAMA-C-2

```
nohic-refpick.sh -s CAMA-C-2-hifi_reads.FOR_ASM.trimmed.fastq.gz -g arabidopsis_pgMC.full.gbz \
                 -i arabidopsis_pgMC.full.hapl -o CAMA-C-2 -t 100 -m 182
```

The `arabidopsis_pgMC.full.gbz` and `arabidopsis_pgMC.full.hapl` files were created using [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md). The pangenome graph contains TAIR10.1 as the reference genome and 47 other assemblies (not including CAMA-C-2) from [Wlodzimierz et al. (2023)](https://doi.org/10.1038/s41586-023-06062-z).

**Outputs**

`CAMA-C-2.synref.fasta`: this is the synref for CAMA-C-2 that will be output to your current working directory.

nohic-refpick.sh will also create a supplemental directory called `CAMA-C-2.refpick_outdir` containing kmers from CAMA-C-2's HiFi reads in [KFF](https://github.com/Kmer-File-Format/kff-reference) format (`CAMA-C-2.kff`) and a personalized graph for CAMA-C-2 (`CAMA-C-2.gbz`).

**Note 1:** If you use the `draft` and `luck` presets of `nohic-asm.sh` (in Subsection 3.3 below), patching of the synref is not necessary. However, if you want to use the other correction presets (e.g. `standard`), we highly recommend patching the synref using a highly contiguous reference genome.    

**Note 2:** If you have large contigs that exceed the limit of the bai index, patching the synref will not work. In that case, please use our [GPatch fork](https://github.com/andyngh/GPatch) 

### 3.3. nohic-asm.sh: Scaffolding Cleaned Contigs Based on a Reference Genome

Once we have clean contigs and synref ready, we can now start with contig error correction and scaffolding using `nohic-asm.sh`. The general usage of the subscript is as following (Please read **our paper** to understand the options for contig correction).

```
nohic-asm.sh -c <contigs.fa> -r <ref.fa> -o <outdir> [options]
Required:
  -c, --contigs <fasta>                               FASTA file containing contig assembly
  -r, --reference <fasta[.gz]>                        FASTA file containing reference genome file
  -o, --output <dir>                                  Common output directory
Optional:
  -fq, --reads <fastq[.gz]>                           FASTQ file with long reads. If not provided, steps 1–2 are auto-disabled; step 3 runs with 'raw' preset.
  -cov, --coverage <int>                              Sequencing coverage for CRAQ (default: 10)
  -t, --threads <int>                                 Thread number (default: 1)
  --ignore-het <yes|no>                               Use sms_clip_coverRate of 0.55 for CRAQ to break heterozygous chimeric contigs (default: no)
  --run-craq <yes|no>                                 Run CRAQ (default: yes)
  --run-inspector <yes|no>                            Run Inspector (default: yes)
  --run-ragtag-correct <yes|no>                       Run RagTag correct (default: yes)
  --run-gap-closing <yes|no>                          Run gap closing with TGSGapCloser (default: no)
  -p, --presets <draft|luck|standard|aggressive|raw>  RagTag correct preset to be run (default: standard)
  --craq-params "<args>"                              User's customized parameters for CRAQ
  --inspector-params "<args>"                         User's customized parameters for inspector.py
  --inspector-correct-params "<args>"                 User's customized parameters for inspector-correct.py
  --ragtag-correct-params "<args>"                    User's customized parameters for RagTag correct
  --ragtag-scf-params "<args>"                        User's customized parameters for RagTag scaffold
  --tgsgapcloser-params "<args>"                      User's customized parameters for TGSGapCloser
  --resume                                            Resume the pipeline at the earliest failed step
  -h, --help                                          Display this help message
```

We will run the following command to correct and scaffold the clean *A. thaliana* CAMA-C-2 contigs.

```
nohic-asm.sh -c CAMA-C-2.asm.bp.p_ctg.pure.for_asm.fa -r CAMA-C-2.synref.fasta -o CAMA-C-2_syn_ref_luck.asm -fq CAMA-C-2-hifi_reads.FOR_ASM.trimmed.fastq.gz \
             -cov 63 -t 100 --run-gap-closing yes -p luck --craq-params "-t 94 -x map-hifi"  --inspector-params "--datatype hifi" \
             --inspector-correct-params "--datatype pacbio-hifi" --ragtag-correct-params "-T corr" --ragtag-scf-params "-C -r -g 2" --tgsgapcloser-params "--tgstype pb"
```

**Note 1:** Please use our **Inspector fork** if you have contigs that are longer than the limit of the bai index

**Note 2:** We recommend users to run CRAQ with 5-6 fewer threads compared to the thread number of the whole subscript. CRAQ run samtools (requiring 5 threads) in parallel with its main command. Thus, setting CRAQ's thread number equal to `nohic-asm` (100 in this example) will crash the subscript.

**Outputs**

`nohic-asm` has 5 steps and the outputs of them are divided into 5 sub-directories, including `1_CRAQ`, `2_Inspector`, `3_RagTag_correct`, `4_Scaffolding`, `5_Gap_closing`. The main outputs of each step are FASTA files in a particular state of contig correction (e.g. `CAMA-C-2.asm.bp.p_ctg.pure.for_asm.craq.inspector.corrected.fa` contains the contigs corrected by CRAQ, Inspector, and RagTag correct). If you do not execute gap closing, you should collect the final scaffolded assembly in `4_Scaffolding`. Otherwise, the final assembly will be in `5_Gap_closing`. All other outputs of CRAQ, Inspector, RagTag correct, RagTag scaffold, and TGSGapcloser are combined together in a directory inside the sub-directory of each step. 

### 3.4. nohic-eval.sh: Assembly Evaluation and Visualization

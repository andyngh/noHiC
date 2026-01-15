# noHiC
Pangenome-based contig scaffolding pipeline with no HiC data 

## [1. Introduction](#1.-Introduction)

### [1.1. noHiC workflow](#1.1.-noHiC-workflow)

noHiC is a reference-guided genome assembly pipeline covering essential steps in contig scaffolding with four independent sub-scripts, including contig decontamination (`nohic-clean`), misassembly correction, reference-guided contig scaffolding (`nohic-refpick` and `nohic-asm`), and assembly evaluation (`nohic-eval`). Given an input FASTA file containing contig-level assembly, the pipeline starts with removing contigs not originating from a specified taxonomic group and optionally organellar contigs. After that, misassemblies in the selected contigs will be corrected in three different approaches using **error-corrected** longs reads (each can be skipped), including clip-based chimeric contig breakings by [CRAQ](https://github.com/JiaoLaboratory/CRAQ/tree/main), small misassembly (base substitutions, expansions, collapses, and haplotype switches) corrections by [Inspector](https://github.com/Maggi-Chen/Inspector), and reference-guided misjoin breakings by [RagTag correct](https://github.com/malonge/RagTag/wiki/correct). The corrected contigs are then scaffolded based on a reference genome with [RagTag scaffold](https://github.com/malonge/RagTag/wiki/scaffold). The most important feature of noHiC is that users can utilized publicly available pangenome graphs built by tools like [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) to synthesized a genetically closed, personalized reference genome (synref) for reference-guided contig error correction and scaffolding based on [KMC](https://github.com/refresh-bio/KMC?tab=readme-ov-file) and [vg haplotype](https://github.com/vgteam/vg/wiki/Haplotype-Sampling), minimizing the amount of false contig correction. Quality of the scaffolded assembly can be assessed based on both metrics (N50, auN, gap number, BUSCO...) and visualizations (i.e. comparing target assembly with a reference genome using a dot plot built by [paf2dotplot](https://github.com/moold/paf2dotplot) and visualizing the positions of remaining misassemblies on chromosomes by `nohic-viz.R`). Please find the detailed descriptions of noHiC in **our paper**.  

![noHiC workflow](https://github.com/andyngh/noHiC/blob/main/Fig1.png)

### [1.2. Dependency citations](#1.2.-Dependency-citations)

Please cite the following tools along with noHiC when you use the pipeline.

**nohic-clean**

[Kraken2](https://doi.org/10.1186/s13059-019-1891-0), [Taxonkit](https://doi.org/10.1016/j.jgg.2021.03.006), [BLAST+](https://doi.org/10.1186/1471-2105-10-421), [Biopython](https://doi.org/10.1093/bioinformatics/btp163)

**nohic-refpick**

[KMC 3](https://doi.org/10.1093/bioinformatics/btx304), [VG](https://doi.org/10.1038/nbt.4227), [VG haplotypes](https://doi.org/10.1038/s41592-024-02407-2), [minimap2](https://doi.org/10.1093/bioinformatics/bty191), [samtools](https://doi.org/10.1093/bioinformatics/btp352), [GPatch](https://doi.org/10.1101/2025.05.22.655567)

**nohic-asm**

[CRAQ](https://doi.org/10.1038/s41467-023-42336-w), [Inspector](https://doi.org/10.1186/s13059-021-02527-4), [RagTag](https://doi.org/10.1186/s13059-022-02823-7), [seqkit](https://doi.org/10.1002/imt2.191), [TGS-GapCloser](https://doi.org/10.1093/gigascience/giaa094)

**nohic-eval**

[gfastats](https://doi.org/10.1093/bioinformatics/btac460), [bioawk](https://github.com/lh3/bioawk), [BUSCO](https://doi.org/10.1093/nar/gkae987), [ggplot2](https://doi.org/10.1007/978-3-319-24277-4), [readr](https://readr.tidyverse.org/), [dplyr](https://dplyr.tidyverse.org/), [paf2dotplot](https://github.com/moold/paf2dotplot?tab=readme-ov-file)

### [1.3. noHiC citation](#1.3.-noHiC-citation)

## [2. Installation](#2.-Installation)

## [3. Tutorial](#3.-Tutorial)

To guide users through the sub-scripts of noHiC, we have provided a small example of *A. thaliana* CAMA-C-2 contig scaffolding, with the essential files in https://github.com/andyngh/noHiC/blob/main/example. In this example, we will try to replicate the manually curated CAMA-C-2 assembly from [this paper](https://doi.org/10.1038/s41586-023-06062-z) using noHiC.

### [3.1. Preparations of HiFi Reads](#3.1.-Preparations-of-HiFi-Reads)

The FASTQ file containing HiFi reads of *A. thaliana* CAMA-C-2 is downloaded from NCBI SRA using `prefetch` and `fastq-dump` from [SRA Toolkit](https://github.com/ncbi/sra-tools). We highly recommend [TGSFilter](https://github.com/HuiyangYu/TGSFilter) for trimming adapters in the reads. 

```
prefetch --max-size 200G ERR10084604

fastq-dump --origfmt ./ERR10084604

tgsfilter -i ERR10084604.fastq -o CAMA-C-2-hifi_reads.ALL.trimmed.fastq.gz -x hifi -t 24
```

In the benchmarks of **our paper**, 90% of the CAMA-C-2 reads will be used for assembly and 10% will be used **only** for assembly evaluation. Though **you do not need to do this step in your own assembly project**, we will conduct the read subsampling here so that the tutorial's results will be the same as in our paper.   

```
# 10% for evaluation
seqkit sample -p 0.10 -j 100 -s 3108 CAMA-C-2-hifi_reads.ALL.trimmed.fastq.gz -o CAMA-C-2-hifi_reads.EVAL.trimmed.fastq.gz

# 90% for assembly
seqkit grep -v -f <(zcat CAMA-C-2-hifi_reads.EVAL.trimmed.fastq.gz | sed -n '1~4s/^@//p') \
  CAMA-C-2-hifi_reads.ALL.trimmed.fastq.gz -o CAMA-C-2-hifi_reads.FOR_ASM.trimmed.fastq.gz
```

### [3.2. nohic-clean.sh: Contaminant Contig Removal](#32-nohic-cleansh-contaminant-contig-removal)

The `nohic-clean` sub-script will be used to remove contaminant contigs from an input contig-level assembly. The general usage of the script is as following:

```
nohic-clean.sh -i <contig.fa> -a <adapter.txt> -o <common_out> -t <threads> -kp <kp_prefix> -tg <taxonomic_group> [-ros <mt.cl.fasta>] [-ioc yes|no] [-m yes|no] [-kdb <kraken2_db_dir>]
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

**Note:** If you have root privilege on your machine and many contig-level assemblies to decontaminate, we recommend executing the following commands before running `nohic-clean` with `-m yes` and **no** `-kdb`.

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

### [3.3. nohic-refpick.sh: Creating Personalized Reference for the Target Genome](#3.3.-nohic-refpick.sh:-Creating-Personalized-Reference-for-the-Target-Genome)

In parallel with contig decontamination, we will create a personalized reference (synref) for *A. thaliana* CAMA-C-2 here using `nohic-refpick.sh`. The general usage of the subscript is as following.

```
nohic-refpick.sh -s <reads.fastq.gz|target_contigs.fa> -g <pangenome_graph.gbz> -i <pangenome_graph.hapl> -o <outprefix> [-t <threads>] [-m <memory>] [-v <0|1|2|3>] [-p <yes|no>] [-r <ref.fa>]
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
                 -i arabidopsis_pgMC.full.hapl -o CAMA-C-2 -t 100 -m 182 -p no
```

The `arabidopsis_pgMC.full.gbz` and `arabidopsis_pgMC.full.hapl` files were created using [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md). The pangenome graph contains TAIR10.1 as the reference genome and 47 other assemblies (not including CAMA-C-2) from [Wlodzimierz et al. (2023)](https://doi.org/10.1038/s41586-023-06062-z).

**Outputs**

`CAMA-C-2.synref.fasta`: this is the synref for CAMA-C-2 that will be output to your current working directory.

nohic-refpick.sh will also create a supplemental directory called `CAMA-C-2.refpick_outdir` containing kmers from CAMA-C-2's HiFi reads in [KFF](https://github.com/Kmer-File-Format/kff-reference) format (`CAMA-C-2.kff`) and a personalized graph for CAMA-C-2 (`CAMA-C-2.gbz`).

**Note 1:** If you use the `draft` and `luck` presets of `nohic-asm.sh` (in Subsection 3.3 below), patching of the synref is not necessary. However, if you want to use the other correction presets (e.g. `standard`), we highly recommend patching the synref using a highly contiguous reference genome.    

**Note 2:** If you have large contigs that exceed the limit of the bai index, patching the synref will not work. In that case, please use our [GPatch fork](https://github.com/andyngh/GPatch) 

### [3.4. nohic-asm.sh: Scaffolding Cleaned Contigs Based on a Reference Genome](#3.4.-nohic-asm.sh:-Scaffolding-Cleaned-Contigs-Based-on-a-Reference-Genome)

Once we have clean contigs and synref ready, we can now start with contig error correction and scaffolding using `nohic-asm.sh`. The general usage of the subscript is as following (Please read **our paper** to understand the options for contig correction).

```
nohic-asm.sh -c <contigs.fa> -r <ref.fa> -o <outdir> -fq <reads.fastq.gz> -cov <sequencing_coverage> -t <threads> [Options]
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

### [3.5. nohic-eval.sh: Assembly Evaluation and Visualization](#3.5.-nohic-eval.sh:-Assembly-Evaluation-and-Visualization)

Quality of the final output assembly of `nohic-asm.sh` (from the `5_Gap_closing` directory) will be evaluated using `nohic-eval.sh`. The general usage of the subscript is as follows:

```
nohic-eval.sh -i <scaffold.fa> -o <outdir> -t <threads> -b <BUSCO_lineage> -r <reads.fastq.gz> --coverage <sequencing_coverage> --sequencing-platform <pb|hifi|ont> --chr-name <chr_names.txt> --reference <ref.fa> [Options]
Required:
      --input-assembly, -i <.fasta|.fa|.fna>    Assembly file (can be gzipped)
      --output-directory, -o <dir>              Common output directory containing all evaluation results
      --threads, -t <num>                       Number of threads

Optional evaluations (default: yes):
      --run-busco <yes|no>                      Run BUSCO
      --run-craq <yes|no>                       Run CRAQ
      --run-inspector <yes|no>                  Run Inspector
      --visualization, -v <yes|no>              Visualize assembly errors

Other optional args:
      --busco-db, -b <str>                      BUSCO lineage database to be used (required if --run-busco yes)
      --reads, -r <fasq[.gz]>                   Long reads for CRAQ/Inspector
      --coverage <num>                          Sequencing coverage for CRAQ
      --sequencing-platform, -p <pb|hifi|ont>   Sequencing platform
      --chr-name <.txt>                         Text file of chromosome names (one per line) for visualization
      --scale-factor <int>                      Stretch misassemblies to this length (bp) for visualization [default: 500000]
      --mm2-params "<str>"                      Minimap2 params for dot plot [default: "-cx asm10"]
      --reference <fasta>                       Reference genome FASTA for dot plot
      --resume                                  Resume nohic-eval from the earliest failed/missing step
      --help, -h                                Show this help
      --version                                 Show version
```

The following command will be run to evaluate the quality of the *A. thaliana* CAMA-C-2 scaffolded assembly.

```
nohic-eval.sh -i CAMA-C-2.asm.bp.p_ctg.pure.for_asm.craq.inspector.corrected.scf.tgs.fa \
              -o CAMA-C-2.asm.bp.p_ctg.pure.for_asm.craq.inspector.corrected.scf.tgs.eval \
              -t 100 -b embryophyta_odb12 -r CAMA-C-2-hifi_reads.EVAL.trimmed.fastq.gz --coverage 7 -p hifi --chr-name ath_chr_names.txt \
              --scale-factor 200000 --reference GCA_946406975.1_CAMA-C-2.PacbioHiFiAssembly_genomic.ed.SELECTED.fa
```

Where, `ath_chr_names.txt` contains chromosome names (one per line) of the `CAMA-C-2.asm.bp.p_ctg.pure.for_asm.craq.inspector.corrected.scf.tgs.fa` assembly (without `>`); `GCA_946406975.1_CAMA-C-2.PacbioHiFiAssembly_genomic.ed.SELECTED.fa` is the public assembly of CAMA-C-2 that we use as the control for structural correctness evaluation.  

**Outputs**

`nohic-eval.sh` stores outputs from each evaluation step as separate sub-directories inside `CAMA-C-2.asm.bp.p_ctg.pure.for_asm.craq.inspector.corrected.scf.tgs.eval`.

*1_Assembly_statistics*

`assembly_stats.txt`: this file contains continuity metrics (N50, total scaffold/contig length, gap number...).

`scaffold_lengths.txt`: this file contains the lengths of scaffolds in the input assembly.

`stats.log`: the log file of this step.

`step_1_done.txt`: step completion marker

*2_BUSCO*

`busco_downloads`: this directory contain the downloaded BUSCO lineage database.

`busco.log`: the log file of this step.

`<input_assembly_name>.busco`: this directory contains the BUSCO results.

`step_2_done.txt`: step completion marker

*3_CRAQ*

`All_CRAQ_outputs`: this directory contains all outputs from CRAQ.

`CRAQ_AQI_metrics.txt`: this file contains the R-AQI and S-AQI metrics showing structural correctness of each scaffold and the whole assembly.

`craq.log`: the log file of this step.

`CSE.csv`: this file contains the CSEs (Clip-based Structural Errors). It has 4 columns, including scaffold name, starting position of a CSE, ending position of a CSE, and error type ("CSE"). This file will be combined with the errors identified by Inspector for misassembly visualization.

`step_3_done.txt`: step completion marker

*4_Inspector*

`All_inspector_outputs`: this directory contains all outputs from Inspector.

`Assembly_statistics_Inspector.txt`: assembly quality metrics calculated by Inspector. You can find the quality value (QV) here.

`small_scale_error.csv`: small scale errors identified by Inspector (BaseSubstitution, SmallCollapse, SmallExpansion...)

`structural_error.csv`: larger structural errors identified by Inspector

`inspector_errors.csv`: this file contains all Inspector-identified assembly errors from `small_scale_error.csv` and `structural_error.csv`. This is the file that will be combined with `3_CRAQ/CSE.csv` for assembly error visualization. It has the same 4 columns as `3_CRAQ/CSE.csv`.

`inspector.log`: the log file of this step.

`step_4_done.txt`: step completion marker

*5_Error_visualization*

`nohic-viz.R`: the helper R script used for visualizing the misassemblies indetified by CRAQ and Inspector on chromosomes. You can use it separately as follows.

```
./nohic-viz.R <chrom_lengths.csv> <misassemblies.csv> <output_prefix> <scale_factor>
```

Where, `chrom_lengths.csv` contains 2 columns **with** headers (`chrom` and `length`). `misassemblies.csv` is the assembly error file with misassemblies from both CRAQ and Inspector. 

`paf2dotplot.R`: this the script used for constructing dot plot showing whole genome alignments between the input assembly and a reference genome. The script usage is in [this repository](https://github.com/moold/paf2dotplot?tab=readme-ov-file). 

`chr_len.csv`: this file contains the chromosome lengths for misassembly visualizations.

`error_to_plots.csv`: this file contains the combined misassemblies indetified by both CRAQ and Inspector.

`plotted_errors.svg`: this file contains the  misassembly visualizations.

`query_to_reference.paf`: the PAF file containing the alignments between the input and the reference.

`query_to_reference.paf.svg`: the dot plot generated from `query_to_reference.paf`.

`step_5_done.txt`: step completion marker

**Evaluation results**

*Continuity metrics*

The continuity metrics of the CAMA-C-2 assembly (after gap closing) are as follows.

```
# scaffolds     6
Total scaffold length   138766131
Average scaffold length 23127688.50
Scaffold N50    26621518
Scaffold auN    27076433.31
Scaffold L50    3
Largest scaffold        33491627
Smallest scaffold       3638613
# contigs       31
Total contig length     138760864
Average contig length   4476156.90
Contig N50      22154653
Contig auN      19037802.18
Contig L50      3
Largest contig  29953688
Smallest contig 105818
# gaps in scaffolds     25
Total gap length in scaffolds   5267
Average gap length in scaffolds 210.68
Gap N50 in scaffolds    679
Gap auN in scaffolds    551.59
Gap L50 in scaffolds    3
Largest gap in scaffolds        1035
Smallest gap in scaffolds       100
Base composition (A:C:G:T)      43769215:25590160:25370347:44031142
GC content %    36.73
# soft-masked bases     0
# segments      31
Total segment length    138760864
Average segment length  4476156.90
# gaps  25
# paths 6
```

In **our paper**, to highlight the benefit of synref in preventing contig fragmentations (in comparison with the assembly guided by a conventional reference - TAIR10.1), we calculated the in-chromosome continuity metrics before gap closing as follows.

```
# scaffolds: 5
Total scaffold length: 135127518
Average scaffold length: 27025503.60
Scaffold N50: 26621518
Scaffold auN: 27707549.45
Scaffold L50: 3
Largest scaffold: 33491627
Smallest scaffold: 21755102
# contigs: 13
Total contig length: 135123951
Average contig length: 10394150.08
Contig N50: 22154653
Contig auN: 19544014.65
Contig L50: 3
Largest contig: 29953688
Smallest contig: 338670
# gaps in scaffolds: 8
Total gap length in scaffolds: 3567
Average gap length in scaffolds: 445.88
Gap N50 in scaffolds: 1013
Gap auN in scaffolds: 766.81
Gap L50 in scaffolds: 2
Largest gap in scaffolds: 1035
Smallest gap in scaffolds: 100
Base composition (A:C:G:T): 42887133:24642914:24453874:43140030
GC content %: 36.33
# soft-masked bases: 0
# segments: 13
Total segment length: 135123951
Average segment length: 10394150.08
# gaps: 8
# paths: 5
```

*Gene space completeness*

The results from BUSCO completeness evaluation of the CAMA-C-2 assembly are as follows.

```
    -------------------------------------------------------------------------------------------
    |Results from dataset embryophyta_odb12                                                    |
    -------------------------------------------------------------------------------------------
    |C:99.1%[S:95.4%,D:3.8%],F:0.2%,M:0.6%,n:2026,E:9.7%                                       |
    |2008    Complete BUSCOs (C)    (of which 195 contain internal stop codons)                |
    |1932    Complete and single-copy BUSCOs (S)                                               |
    |76    Complete and duplicated BUSCOs (D)                                                  |
    |5    Fragmented BUSCOs (F)                                                                |
    |13    Missing BUSCOs (M)                                                                  |
    |2026    Total BUSCO groups searched                                                       |
    -------------------------------------------------------------------------------------------

```

*Structural correctness*

The AQIs and QV showing the strutural completeness of CAMA-C-2 are as follows.

```
# From CRAQ:

Short Report:
#Chr    Covered.Rate    Low-confident.Rate      Avg.CRH Avg.CSH Avg.CRE(R-AQI)  Avg.CSE(S-AQI)
Genome  0.942858357850962       0.00782285262516361     0       0       0.266635654297737(97.3688768987995)     0(100)
chr5_RagTag     0.970727179188941       0.00869195109824369     0       0       0.0687832156297752(99.3145279946857)    0(100)
chr4_RagTag     0.942827738501308       0.0208867110989366      0       0       0.0910201686130419(99.0939281100937)    0(100)
chr1_RagTag     0.973810440811941       0.00536644560843683     0       0       0.122644843406757(98.7810417922147)     0(100)
chr3_RagTag     0.983733272931072       0.00203602206170129     0       0       0.152738943371082(98.4842159968631)     0(100)
chr2_RagTag     0.964547262104562       0.00347490885241959     0       0       0.28593450975051(97.1811473172815)      0(100)
Chr0_RagTag     2.74829845971613e-07    0.0079381874526949      0       0       19000000(0)     0(100)

# From Inspector:

QV      64.83304855386444
```

The following figure shows the misassemblies in CAMA-C-2. 

![Arabidopsis_assembly_errors](https://github.com/andyngh/noHiC/blob/main/figures/error_synref.png)

The resulting CAMA-C-2 assembly (y-axis) has strong strutural agreement to the public, manually curated assembly (x-axis).

![Arabidopsis_dot_plot](https://github.com/andyngh/noHiC/blob/main/figures/ara_synref_to_NCBI.paf.png)




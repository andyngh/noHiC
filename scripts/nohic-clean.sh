#!/usr/bin/env bash
set -euo pipefail
#-------------------#
# nohic-clean v1.0  #
#-------------------#

usage() {
  cat <<EOF
nohic-clean v1.0 - Contig decontamination pipeline
___________________________________________________

Usage: $0 -i <contig.fa> -a <adapter.txt> -o <common_out> -t <threads> -kp <kp_prefix> -tg <taxonomic_group> [-ros <mt.cl.fasta>] [-ioc yes|no] [-m yes|no] [-kdb <kraken2_db_dir>] [--resume]
Required:
  -i, --input-fasta <.fa|.fasta>                                FASTA file with raw input contigs 
  -a, --adapter-sequence <.txt>               	                A text file containing adapter sequences (one adapter per line)
  -o, --output-directory <dir>               	                Common output directory for storing results from contig cleaning steps
  -t, --threads <int>                        	                Thread number for kraken2 and blastn
  -kp, --prefix-for-taxonomic-classification <str>  	        Prefix for kraken2 outputs (e.g. sample_1)
  -tg, --taxonomic-group <str>               		        Include contigs from this taxonomic group (e.g. Viridiplantae, Bacteria, Mammalia...)
Optional:
  -ros, --reference-organellar-sequences <.fa|.fasta>  	        FASTA file containing reference organellar DNA sequences 
  -ioc, --identify-organellar-contigs <yes|no>    	        Check if there is organellar DNA in the contigs (Default: no)
  -m,  --memory-mapping <yes|no>                                Use --memory-mapping option of kraken2 (Default: yes)
  -kdb, --kraken2-database <dir>                                Path to the directory containing a kraken2 database (*.k2d files) (required if -m no)
  --resume                            		                Resume the pipeline from the earliest failed steps
  -h|--help                           		                Display this help message
EOF
  exit 1
}

# Defaults
IDENTIFY_ORGANELLAR="no"
RESUME="no"
REFERENCE_ORG=""
THREADS=1
MEMORY_MAPPING="yes"
KRAKEN_DB=""

# Determine script directory
SCRIPT_DIR="$(dirname "$(realpath "$0")")"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-fasta) INPUT_FASTA="$2"; shift 2;;
    -a|--adapter-sequence) ADAPTER="$2"; shift 2;;
    -o|--output-directory) OUTDIR="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    -kp|--prefix-for-taxonomic-classification) KP="$2"; shift 2;;
    -tg|--taxonomic-group) TAXGROUP="$2"; shift 2;;
    -ros|--reference-organellar-sequences) REFERENCE_ORG="$2"; shift 2;;
    -ioc|--identify-organellar-contigs) IDENTIFY_ORGANELLAR="$2"; shift 2;;
    -m|--memory-mapping) MEMORY_MAPPING="$2"; shift 2;;
    -kdb|--kraken2-database) KRAKEN_DB="$2"; shift 2;;
    --resume) RESUME="yes"; shift;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

# Required checks
: "${INPUT_FASTA:?Missing -i/--input-fasta}"
: "${ADAPTER:?Missing -a/--adapter-sequence}"
: "${OUTDIR:?Missing -o/--output-directory}"
: "${THREADS:?Missing -t/--threads}"
: "${KP:?Missing -kp/--prefix-for-taxonomic-classification}"
: "${TAXGROUP:?Missing -tg/--taxonomic-group}"

# Validate memory mapping options and conditional requirement
if [[ "$MEMORY_MAPPING" != "yes" && "$MEMORY_MAPPING" != "no" ]]; then
  echo "Error: -m/--memory-mapping must be 'yes' or 'no' (got: $MEMORY_MAPPING)" >&2
  exit 1
fi
if [[ "$MEMORY_MAPPING" == "no" ]]; then
  : "${KRAKEN_DB:?Missing -kdb/--kraken2-database when -m no}"
  if [[ ! -d "$KRAKEN_DB" ]]; then
    echo "Error: Provided kraken2 database directory does not exist: $KRAKEN_DB" >&2
    exit 1
  fi
fi

# Normalize paths
INPUT_FASTA=$(realpath "$INPUT_FASTA")
ADAPTER=$(realpath "$ADAPTER")
OUTDIR=$(realpath "$OUTDIR")
THREADS=${THREADS}
KP=${KP}
TAXGROUP=${TAXGROUP}
if [[ -n "$REFERENCE_ORG" ]]; then REFERENCE_ORG=$(realpath "$REFERENCE_ORG"); fi
if [[ -n "$KRAKEN_DB" ]]; then KRAKEN_DB=$(realpath "$KRAKEN_DB"); fi

# Derived names & step dirs
BASENAME=$(basename "$INPUT_FASTA")
BASENAME_NO_EXT="${BASENAME%.*}"

STEP1_DIR="$OUTDIR/1_adapter_content_check"
STEP2_DIR="$OUTDIR/2_contig_taxonomic_classification"
STEP3_DIR="$OUTDIR/3_organellar_contigs_identification"
STEP4_DIR="$OUTDIR/4_contig_purification"

mkdir -p "$STEP1_DIR" "$STEP2_DIR" "$STEP3_DIR" "$STEP4_DIR"

# Helper for step done files
done_file() {
  local stepdir="$1"
  local stepnum="$2"
  echo "$stepdir/step_${stepnum}_done.txt"
}

# Locations of helper scripts in script dir (if present)
FILTER_SCRIPT_IN_SCRIPT_DIR="$SCRIPT_DIR/filter_blast.py"
REMOVE_SCRIPT_IN_SCRIPT_DIR="$SCRIPT_DIR/remove_contamination.py"

####################################################################################################
# Resume logic:								                           #
# - Determine earliest step whose "step_*_done.txt" is missing OR does not contain "ok" to restart.#
# - START_STEP=5 means all steps 1..4 are ok.						           #
####################################################################################################
START_STEP=1
if [[ "$RESUME" == "yes" ]]; then
  declare -A __STEP_DIRS
  __STEP_DIRS[1]="$STEP1_DIR"
  __STEP_DIRS[2]="$STEP2_DIR"
  __STEP_DIRS[3]="$STEP3_DIR"
  __STEP_DIRS[4]="$STEP4_DIR"

  START_STEP=5
  for s in 1 2 3 4; do
    df="${__STEP_DIRS[$s]}/step_${s}_done.txt"
    if [[ ! -f "$df" ]]; then
      START_STEP=$s
      break
    fi
    
    status=$(tr -d '\n\r\t ' < "$df" 2>/dev/null | tr '[:upper:]' '[:lower:]')
    if [[ "$status" != "ok" ]]; then
      START_STEP=$s
      break
    fi
  done

  if (( START_STEP < 5 )); then
    echo "Resume requested: nohic-clean will restart from step: $START_STEP."
  else
    echo "Resume requested: all steps appear complete (1–4). Nothing to re-run."
  fi
fi


################################################################################
# Step 1: Adapter content check                                                #
################################################################################
STEP1_DONE=$(done_file "$STEP1_DIR" 1)
STEP1_LOG="$STEP1_DIR/adapter_check.log"
if [[ "$RESUME" == "yes" && $START_STEP -gt 1 ]]; then
  echo "Step 1 already done. Skipping step 1."
else
  trap 'echo "failed" > "$STEP1_DONE"' ERR
  trap 'echo "failed" > "$STEP1_DONE"; exit 130' INT TERM
  echo "Running adapter check..." | tee "$STEP1_LOG"
  if [[ ! -f "$ADAPTER" ]]; then
    echo "Adapter file not found: $ADAPTER" | tee -a "$STEP1_LOG" >&2
    echo "failed" > "$STEP1_DONE"

    exit 1
  fi
  if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "Input fasta not found: $INPUT_FASTA" | tee -a "$STEP1_LOG" >&2
    echo "failed" > "$STEP1_DONE"

    exit 1
  fi

  echo "Searching for adapter occurrences..." | tee -a "$STEP1_LOG"
  FGREP_OUT="$STEP1_DIR/fgrep_matches.txt"
  : > "$FGREP_OUT"
  # run fgrep with line numbers; store to a file (filename:lineno:content)
  fgrep -n -H -f "$ADAPTER" "$INPUT_FASTA" > "$FGREP_OUT" 2>>"$STEP1_LOG" || true

  if [[ ! -s "$FGREP_OUT" ]]; then
    echo "Adapter check: OK" | tee -a "$STEP1_LOG"
    echo "None" > "$STEP1_DIR/adapter_content.txt"
    echo "ok" > "$STEP1_DONE"
  else
    echo "Adapter matches detected. Preparing $STEP1_DIR/adapter_matches.fa" | tee -a "$STEP1_LOG"
    OUT_FA="$STEP1_DIR/adapter_matches.fa"
    : > "$OUT_FA"
    OUT_NAMES="$STEP1_DIR/adapter_matches.txt"
    : > "$OUT_NAMES"

    # collect unique matched line numbers (2nd field)
    mapfile -t linenums < <(cut -d: -f2 "$FGREP_OUT" | sort -n | uniq)

    declare -A seen_pairs

    for ln in "${linenums[@]}"; do
      # ensure ln is a positive integer
      if ! [[ "$ln" =~ ^[0-9]+$ ]]; then
        continue
      fi
      seqline_num=$((ln))
      # search backward from seqline_num-1 until find a header starting with '>'
      probe=$((seqline_num - 1))
      header_line=""
      while [[ $probe -ge 1 ]]; do
        h=$(sed -n "${probe}p" "$INPUT_FASTA" 2>/dev/null || true)
        # skip empty lines
        if [[ -z "$h" ]]; then
          probe=$((probe - 1))
          continue
        fi
        if [[ "${h:0:1}" == ">" ]]; then
          header_line="$h"
          break
        fi
        probe=$((probe - 1))
      done

      # if didn't find a header, warn and skip
      if [[ -z "$header_line" ]]; then
        echo "Warning: could not find header for matched line $seqline_num (skipping)" | tee -a "$STEP1_LOG"
        continue
      fi

      seq_line_content=$(sed -n "${seqline_num}p" "$INPUT_FASTA" 2>/dev/null || true)

      # dedupe based on header + seq line
      key="${header_line}||${seq_line_content}"
      if [[ -z "${seen_pairs[$key]:-}" ]]; then
        printf "%s\n" "$header_line" >> "$OUT_FA"
        printf "%s\n" "$seq_line_content" >> "$OUT_FA"
        # write first token of header to names file
        first_token=$(printf "%s" "$header_line" | awk '{print $1}')
        printf "%s\n" "$first_token" >> "$OUT_NAMES"
        seen_pairs["$key"]=1
      fi
    done

    MATCHES=$(grep -c '^>' "$OUT_FA" || true)
    echo "Adapter matches found in contigs: $MATCHES" | tee -a "$STEP1_LOG"
    echo "Adapter check: failed. Adapter sequence detected in contigs. See $OUT_FA and $STEP1_LOG" | tee -a "$STEP1_LOG" >&2
    echo "failed" > "$STEP1_DONE"

    exit 2
  fi
fi

################################################################################
# Step 2: Contig taxonomic classification with kraken2 + taxonkit              #
################################################################################
STEP2_DONE=$(done_file "$STEP2_DIR" 2)
KR_OUTPUT="$STEP2_DIR/${KP}.kr"
REPORT_FILE="$STEP2_DIR/${KP}.report"
LINEAGE_FILE="$STEP2_DIR/${KP}.lineage"
CONTAM_FILE="$STEP2_DIR/name_of_contaminant_contigs.txt"
STEP2_LOG="$STEP2_DIR/kraken2_and_taxonkit.log"

if [[ "$RESUME" == "yes" && $START_STEP -gt 2 ]]; then
  echo "Step 2 already done. Skipping step 2."
else
  trap 'echo "failed" > "$STEP2_DONE"' ERR
  trap 'echo "failed" > "$STEP2_DONE"; exit 130' INT TERM

  if [[ "$MEMORY_MAPPING" == "yes" ]]; then
    echo "Checking for kraken2 database files (*.k2d) in /dev/shm/ ..." | tee "$STEP2_LOG"
    K2D_COUNT=$(sh -c 'ls /dev/shm/*.k2d 2>/dev/null | wc -l' || true)
    if [[ "$K2D_COUNT" -eq 0 ]]; then
      cat <<EOF | tee -a "$STEP2_LOG"
Check kraken2 database: Failed. Cannot find *.k2d files in /dev/shm/. Please do the following steps before re-running the pipeline:

If your /dev/shm directory size is smaller than the downloaded kraken2 database, resize it first by this command: 

sudo mount -o remount,size=your_database_size /dev/shm

For example, if your kraken2 database size is 242G, run this:

sudo mount -o remount,size=242G /dev/shm

Copy your *.k2d files to /dev/shm/

cp /path/to/your/database/directory/*.k2d /dev/shm/
EOF
      echo "failed" > "$STEP2_DONE"

      exit 3
    else
      echo "Check kraken2 database: OK (found $K2D_COUNT .k2d files)." | tee -a "$STEP2_LOG"
    fi

    echo "Running kraken2 with --memory-mapping option..." | tee -a "$STEP2_LOG"
    kraken2 --db /dev/shm/ --threads "$THREADS" --memory-mapping --output "$KR_OUTPUT" --report "$REPORT_FILE" "$INPUT_FASTA" >> "$STEP2_LOG" 2>&1
  else
    echo "Running kraken2 using provided database..." | tee "$STEP2_LOG"
    kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --output "$KR_OUTPUT" --report "$REPORT_FILE" "$INPUT_FASTA" >> "$STEP2_LOG" 2>&1
  fi

  echo "kraken2 finished. Outputs: $KR_OUTPUT and $REPORT_FILE" | tee -a "$STEP2_LOG"

  # Taxonkit: check for taxdump in $HOME/.taxonkit
  TAXONKIT_DIR="$HOME/.taxonkit"
  DMP_COUNT=$(sh -c "ls $TAXONKIT_DIR/*.dmp 2>/dev/null | wc -l" || true)
  if [[ "$DMP_COUNT" -eq 0 ]]; then
    echo "Downloading database for Taxonkit..." | tee -a "$STEP2_LOG"
    wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O "$STEP2_DIR/taxdump.tar.gz" >> "$STEP2_LOG" 2>&1
    tar -zxvf "$STEP2_DIR/taxdump.tar.gz" -C "$STEP2_DIR" >> "$STEP2_LOG" 2>&1
    mkdir -p "$TAXONKIT_DIR"
    if [[ -d "$STEP2_DIR/taxdump" ]]; then
      mv "$STEP2_DIR/taxdump/"*.dmp "$TAXONKIT_DIR/" 2>> "$STEP2_LOG" || true
    else
      mv "$STEP2_DIR/"*.dmp "$TAXONKIT_DIR/" 2>> "$STEP2_LOG" || true
    fi
    echo "Taxonkit taxdump files placed in $TAXONKIT_DIR" | tee -a "$STEP2_LOG"
  else
    echo "Taxonkit database already present (found $DMP_COUNT .dmp files)." | tee -a "$STEP2_LOG"
  fi

  echo "Running taxonkit lineage..." | tee -a "$STEP2_LOG"
  cut -f3 "$KR_OUTPUT" | taxonkit lineage > "$LINEAGE_FILE" 2>> "$STEP2_LOG" || true
  echo "Taxonkit lineage written to: $LINEAGE_FILE" | tee -a "$STEP2_LOG"

  echo "Generating list of contaminant contigs..." | tee -a "$STEP2_LOG"
  paste "$KR_OUTPUT" "$LINEAGE_FILE" | cut -f 1,2,3,4,6,7 | grep -v "$TAXGROUP" | cut -f2 | sort | uniq > "$CONTAM_FILE" 2>> "$STEP2_LOG" || true
  echo "Contaminant contigs list: $CONTAM_FILE" | tee -a "$STEP2_LOG"

  echo "Step 2 done! Files written to: $STEP2_DIR" | tee -a "$STEP2_LOG"
  echo "ok" > "$STEP2_DONE"
fi

################################################################################
# Step 3: Identification of organellar contigs (BLAST check)                   #
################################################################################
STEP3_DONE=$(done_file "$STEP3_DIR" 3)
BLAST_OUT="$STEP3_DIR/blast_out.tab"
ORG_CONTS="$STEP3_DIR/organellar_dna_containing_contigs.txt"
STEP3_LOG="$STEP3_DIR/blast_formatdb_and_filter.log"

if [[ "$RESUME" == "yes" && $START_STEP -gt 3 ]]; then
  echo "Step 3 already done. Skipping step 3."
else
  trap 'echo "failed" > "$STEP3_DONE"' ERR
  trap 'echo "failed" > "$STEP3_DONE"; exit 130' INT TERM
  if [[ "$IDENTIFY_ORGANELLAR" == "no" ]]; then
    echo "Organellar contig identification disabled (-ioc no). Creating empty file." | tee "$STEP3_LOG"
    : > "$ORG_CONTS"
    echo "ok" > "$STEP3_DONE"
  else
    if [[ -z "$REFERENCE_ORG" ]]; then
      echo "Error: identify-organellar-contigs=yes but no -ros/--reference-organellar-sequences provided." | tee -a "$STEP3_LOG" >&2
      echo "failed" > "$STEP3_DONE"

      exit 4
    fi
    REFERENCE_ORG=$(realpath "$REFERENCE_ORG")

    echo "Checking for existing BLAST database files..." | tee -a "$STEP3_LOG"
    DB_NHR="${REFERENCE_ORG}.nhr"
    DB_NIN="${REFERENCE_ORG}.nin"
    DB_NSQ="${REFERENCE_ORG}.nsq"
    if [[ -f "$DB_NHR" && -f "$DB_NIN" && -f "$DB_NSQ" ]]; then
      echo "All three BLAST database files (.nhr .nin .nsq) exist for $REFERENCE_ORG. Skipping makeblastdb." | tee -a "$STEP3_LOG"
    else
      echo "Creating BLAST database with makeblastdb..." | tee -a "$STEP3_LOG"
      makeblastdb -in "$REFERENCE_ORG" -dbtype nucl -parse_seqids -out "$REFERENCE_ORG" >> "$STEP3_LOG" 2>&1
      echo "makeblastdb finished" | tee -a "$STEP3_LOG"
    fi

    echo "Running blastn against organellar reference..." | tee -a "$STEP3_LOG"
    blastn -query "$INPUT_FASTA" -db "$REFERENCE_ORG" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -num_threads "$THREADS" -out "$BLAST_OUT" >> "$STEP3_LOG" 2>&1
    echo "BLAST finished. Output: $BLAST_OUT" | tee -a "$STEP3_LOG"

    # Decide where filter_blast.py comes from
    if [[ -f "$FILTER_SCRIPT_IN_SCRIPT_DIR" ]]; then
      FILTER_SCRIPT="$FILTER_SCRIPT_IN_SCRIPT_DIR"
      echo "Using existing filter_blast.py from script directory: $FILTER_SCRIPT" | tee -a "$STEP3_LOG"
    else
      FILTER_SCRIPT="$STEP3_DIR/filter_blast.py"
      if [[ ! -f "$FILTER_SCRIPT" ]]; then
        cat > "$FILTER_SCRIPT" <<'PY'
import sys

for l in open(sys.argv[1]):
    l_arr=l.rstrip().split("\t")
    if(float(l_arr[2])>=90.0 and float(l_arr[12])>=90.0):
        print(l.rstrip())
PY
        chmod +x "$FILTER_SCRIPT"
        echo "Created filter_blast.py at $FILTER_SCRIPT" | tee -a "$STEP3_LOG"
      else
        echo "filter_blast.py already exists at $FILTER_SCRIPT" | tee -a "$STEP3_LOG"
      fi
    fi

    echo "Filtering BLAST results with $FILTER_SCRIPT..." | tee -a "$STEP3_LOG"
    python3 "$FILTER_SCRIPT" "$BLAST_OUT" > "$STEP3_DIR/filter_blast_out.tab" 2>> "$STEP3_LOG" || true
    cut -f 1 "$STEP3_DIR/filter_blast_out.tab" | sort | uniq > "$ORG_CONTS"
    echo "Organellar contigs written to: $ORG_CONTS" | tee -a "$STEP3_LOG"

    echo "Step 3 done! Organellar contigs list at: $ORG_CONTS" | tee -a "$STEP3_LOG"
    echo "ok" > "$STEP3_DONE"
  fi
fi

################################################################################
# Step 4: Removal of contaminant contigs                                       #
################################################################################
STEP4_DONE=$(done_file "$STEP4_DIR" 4)
COMBINED_LIST="$STEP4_DIR/contigs_from_contaminant_n_organelle.txt"
OUTPUT_PURE_FASTA="$STEP4_DIR/${BASENAME_NO_EXT}.pure.fa"
STEP4_LOG="$STEP4_DIR/remove_contamination.log"

if [[ "$RESUME" == "yes" && $START_STEP -gt 4 ]]; then
  echo "Step 4 already done (found $STEP4_DONE). Skipping step 4."
else
  trap 'echo "failed" > "$STEP4_DONE"' ERR
  trap 'echo "failed" > "$STEP4_DONE"; exit 130' INT TERM
  echo "Combining contaminant and organellar lists..." | tee "$STEP4_LOG"
  touch "$CONTAM_FILE" "$ORG_CONTS"
  cat "$CONTAM_FILE" "$ORG_CONTS" | sort | uniq > "$COMBINED_LIST"
  echo "Combined list at: $COMBINED_LIST" | tee -a "$STEP4_LOG"

  # Decide where remove_contamination.py comes from
  if [[ -f "$REMOVE_SCRIPT_IN_SCRIPT_DIR" ]]; then
    REMOVE_SCRIPT="$REMOVE_SCRIPT_IN_SCRIPT_DIR"
    echo "Using existing remove_contamination.py from script directory: $REMOVE_SCRIPT" | tee -a "$STEP4_LOG"
  else
    REMOVE_SCRIPT="$STEP4_DIR/remove_contamination.py"
    if [[ ! -f "$REMOVE_SCRIPT" ]]; then
      cat > "$REMOVE_SCRIPT" <<'PY'
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
c_s=set()
for l in open(sys.argv[1]):
    c_s.add(l.rstrip())
with open(sys.argv[2]) as handle:
    for seq in SeqIO.parse(handle, "fasta"):
        if(seq.id not in c_s):
            sys.stdout.write(seq.format("fasta"))
            sys.stdout.flush()
handle.close()
PY
      chmod +x "$REMOVE_SCRIPT"
      echo "Created remove_contamination.py at $REMOVE_SCRIPT" | tee -a "$STEP4_LOG"
    else
      echo "remove_contamination.py already exists at $REMOVE_SCRIPT" | tee -a "$STEP4_LOG"
    fi
  fi

  echo "Removing contaminant contigs using $REMOVE_SCRIPT..." | tee -a "$STEP4_LOG"
  python3 "$REMOVE_SCRIPT" "$COMBINED_LIST" "$INPUT_FASTA" > "$OUTPUT_PURE_FASTA" 2>> "$STEP4_LOG" || true
  echo "Pure contigs fasta written to: $OUTPUT_PURE_FASTA" | tee -a "$STEP4_LOG"
  echo "ok" > "$STEP4_DONE"
fi

echo "Pipeline finished. Outputs collected under: $OUTDIR"
echo "Step summary:"
for d in "$STEP1_DIR" "$STEP2_DIR" "$STEP3_DIR" "$STEP4_DIR"; do
  echo " - $d"
done

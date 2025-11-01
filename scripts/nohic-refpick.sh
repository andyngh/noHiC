#!/usr/bin/env bash
set -euo pipefail

VERSION="nohic-refpick v1.0"

# Default parameters

THREADS=1
RAM_GB=32
K=29
OUTPREFIX=""
CONTIGS=""
GBZ_IN=""
HAP_INDEX=""
VERBOSITY=2

print_help() {
  cat <<'EOF'
nohic-refpick v1.0 — Build a synthetic reference from a pangenome graph
_______________________________________________________________________

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
  -h, --help                            Display this help message
      --version                         Display version number

EOF
}

# Args
PARSED=$(getopt -o s:g:i:o:t:m:k:v:h \
  --long input-sequence:,gbz:,hapidx:,outprefix:,threads:,ram-gb:,kmer:,verbosity:,help,version \
  -n nohic-refpick -- "$@") || { echo "Error: failed to parse options"; exit 2; }
eval set -- "$PARSED"

while true; do
  case "$1" in
    -s|--input-sequence) CONTIGS=$2; shift 2 ;;
    -g|--gbz)       GBZ_IN=$2; shift 2 ;;
    -i|--hapidx)    HAP_INDEX=$2; shift 2 ;;
    -o|--outprefix) OUTPREFIX=$2; shift 2 ;;
    -t|--threads)   THREADS=$2; shift 2 ;;
    -m|--ram-gb)    RAM_GB=$2; shift 2 ;;
    -k|--kmer)      K=$2; shift 2 ;;
    -v|--verbosity) VERBOSITY=$2; shift 2 ;;
    -h|--help)      print_help; exit 0 ;;
       --version)
         echo "${VERSION}"
         exit 0 ;;
    --) shift; break ;;
    *) echo "Internal parsing error"; exit 2 ;;
  esac
done

# validations
[[ -z "$CONTIGS"   ]] && { echo "Error: --input-sequence is required"; exit 2; }
[[ -z "$GBZ_IN"    ]] && { echo "Error: --gbz is required"; exit 2; }
[[ -z "$HAP_INDEX" ]] && { echo "Error: --hapidx is required"; exit 2; }
[[ -z "$OUTPREFIX" ]] && { echo "Error: --outprefix is required"; exit 2; }

[[ -f "$CONTIGS"   ]] || { echo "Error: input sequence file not found: $CONTIGS"; exit 1; }
[[ -f "$GBZ_IN"    ]] || { echo "Error: GBZ file not found: $GBZ_IN"; exit 1; }
[[ -f "$HAP_INDEX" ]] || { echo "Error: hap index not found: $HAP_INDEX"; exit 1; }

command -v kmc >/dev/null || { echo "Error: kmc not found in PATH"; exit 127; }
command -v vg  >/dev/null || { echo "Error: vg not found in PATH"; exit 127; }

OUTDIR="${OUTPREFIX}.refpick_outdir"
mkdir -p "$OUTDIR"

# infer input mode for KMC
infer_kmc_mode() {
  local f="$1"
  local lf="${f,,}"
  # extension-based guess
  case "$lf" in
    *.fa|*.fna|*.fasta|*.fa.gz|*.fna.gz|*.fasta.gz|*.fa.bz2|*.fna.bz2|*.fasta.bz2|*.fa.xz|*.fna.xz|*.fasta.xz)
      echo "-fm"; return ;;
    *.fq|*.fastq|*.fq.gz|*.fastq.gz|*.fq.bz2|*.fastq.bz2|*.fq.xz|*.fastq.xz)
      echo "-fq"; return ;;
  esac
  # fallback: sniff first byte
  local decomp="cat"
  [[ "$lf" == *.gz  ]] && decomp="gzip -dc"
  [[ "$lf" == *.bz2 ]] && decomp="bzip2 -dc"
  [[ "$lf" == *.xz  ]] && decomp="xz -dc"
  local first_char
  if first_char=$($decomp "$f" 2>/dev/null | head -c 1); then
    if [[ "$first_char" == ">" ]]; then echo "-fm"; return; fi
    if [[ "$first_char" == "@" ]]; then echo "-fq"; return; fi
  fi
  # default to FASTA mode if unknown
  echo "-fm"
}

KMC_MODE=$(infer_kmc_mode "$CONTIGS")
echo "[INFO] Detected input: $CONTIGS.  KMC will be run with mode: ${KMC_MODE#-}"

echo "[INFO] Running KMC: k=${K}, RAM=${RAM_GB}GB, threads=${THREADS}"

kmc -k"$K" -m"$RAM_GB" -okff -t"$THREADS" -hp ${KMC_MODE} \
  "$CONTIGS" "$OUTDIR/$OUTPREFIX" "$OUTDIR"

# vg haplotypes

GBZ_OUT="$OUTDIR/$OUTPREFIX.gbz"
KFF="$OUTDIR/$OUTPREFIX.kff"

echo "[INFO] Running vg haplotypes (verbosity=$VERBOSITY)"
vg haplotypes -v "$VERBOSITY" -t "$THREADS" \
  --preset default --num-haplotypes 1 \
  -i "$HAP_INDEX" -k "$KFF" -g "$GBZ_OUT" "$GBZ_IN"

# extract synthetic reference
FASTA_OUT="${OUTPREFIX}.synref.fasta"
echo "[INFO] Extracting synthetic reference"

vg paths -t "$THREADS" --extract-fasta -x "$GBZ_OUT" --paths-by recombination > "$FASTA_OUT"

echo "[INFO] nohic-refpick completed."
echo "Outputs:"
echo "  $FASTA_OUT"
echo "  $OUTDIR/ (intermediates: .kff, .gbz)"

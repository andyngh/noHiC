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

PATCH="yes"
REFERENCE=""

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
  -p, --patch <yes|no>                  Run personalized reference patching step (default: yes)
  -r, --reference <.fasta>              Reference fasta for patching (required if --patch yes)
  -h, --help                            Display this help message
      --version                         Display version number

EOF
}

# Args
PARSED=$(getopt -o s:g:i:o:t:m:k:v:p:r:h \
  --long input-sequence:,gbz:,hapidx:,outprefix:,threads:,ram-gb:,kmer:,verbosity:,patch:,reference:,help,version \
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
    -p|--patch)     PATCH=$2; shift 2 ;;
    -r|--reference) REFERENCE=$2; shift 2 ;;
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

if [[ "$PATCH" != "yes" && "$PATCH" != "no" ]]; then
  echo "Error: --patch must be \"yes\" or \"no\" (got: $PATCH)"
  exit 2
fi

if [[ "$PATCH" == "yes" && -z "${REFERENCE:-}" ]]; then
  echo "Error: --reference is required when --patch yes"
  exit 2
fi

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

echo "[INFO] personalized reference generation completed."
echo "Outputs:"
echo "  $FASTA_OUT"
echo "  $OUTDIR/ (intermediates: .kff, .gbz)"

if [[ "$PATCH" == "yes" ]]; then

  ASM_DECOMP_SCRIPT="$OUTDIR/Asm_Decomposing.sh"
  cat <<'EOF_ASM_DECOMP' > "$ASM_DECOMP_SCRIPT"
#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <assembly.fasta> [output.fasta]" >&2
    exit 1
fi

IN="$1"
OUT="${2:-/dev/stdout}"

awk -v OFS="" '
BEGIN {
    # No contig yet
}

# Header line: start a new scaffold
/^>/ {
    # Flush any previous scaffold sequence
    if (seq != "") {
        process_seq(scaf_id, seq)
    }
    # Get scaffold ID (up to first whitespace)
    scaf_id = $1
    sub(/^>/, "", scaf_id)
    sub(/ .*/, "", scaf_id)

    # Reset sequence and contig counter for this scaffold
    seq = ""
    ctg_idx = 0
    next
}

# Sequence lines: append (remove whitespace and make uppercase)
{
    gsub(/[ \t\r\n]/, "", $0)
    if ($0 != "") {
        seq = seq toupper($0)
    }
    next
}

END {
    # Flush last scaffold
    if (seq != "") {
        process_seq(scaf_id, seq)
    }
}

# Function to process one scaffold sequence: split on Ns, output contigs
function process_seq(id, s,   i, len, start, in_block, base, seg) {
    len = length(s)
    start = 1
    in_block = 0

    for (i = 1; i <= len; i++) {
        base = substr(s, i, 1)
        if (base != "N") {
            if (!in_block) {
                # Start new contig
                start = i
                in_block = 1
            }
        } else {
            # base == N: end of a contig block if we were in one
            if (in_block) {
                seg = substr(s, start, i - start)
                output_contig(id, seg)
                in_block = 0
            }
        }
    }

    # Trailing contig block to end of sequence
    if (in_block) {
        seg = substr(s, start, len - start + 1)
        output_contig(id, seg)
    }
}

# Output a contig with proper name, skipping empty segments
function output_contig(id, seq_seg,   header) {
    if (length(seq_seg) == 0) {
        return
    }
    ctg_idx++
    header = ">ctg" ctg_idx "_" id
    print header
    # Wrap sequence at 60 bp per line (FASTA formatting)
    line_len = 60
    for (pos = 1; pos <= length(seq_seg); pos += line_len) {
        print substr(seq_seg, pos, line_len)
    }
}
' "$IN" > "$OUT"
EOF_ASM_DECOMP
  chmod +x "$ASM_DECOMP_SCRIPT"

  SYNREF_CTGS="${OUTPREFIX}.synref.ctgs.fasta"
  echo "[INFO] Decomposing personalized reference into contigs: $SYNREF_CTGS"
  "$ASM_DECOMP_SCRIPT" "$FASTA_OUT" "$SYNREF_CTGS"

  echo "[INFO] Patching in progress"
  minimap2 -t "$THREADS" -a -x asm10 "$REFERENCE" "$SYNREF_CTGS" > synref.sam
  samtools view -bS -@ "$THREADS" synref.sam > synref_to_donor_aln.bam
  GPatch -q synref_to_donor_aln.bam -r "$REFERENCE" -x "${OUTPREFIX}.synref"
  rm synref.sam

  echo "[INFO] Patching completed."
fi

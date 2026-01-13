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

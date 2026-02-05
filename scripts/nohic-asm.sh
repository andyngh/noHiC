#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------#
#                       nohic-asm v1.0.0                            #
# ------------------------------------------------------------------#

# --- Traps --------------------------------------------------------
# mark "running" step as "failed" on any exit or kill
mark_running_failed() {
  if [[ -n "${CURRENT_STEP_FILE:-}" && -f "$CURRENT_STEP_FILE" ]]; then
    local s; s="$(tr -d ' \t\r\n' < "$CURRENT_STEP_FILE" | tr '[:upper:]' '[:lower:]')"
    if [[ "$s" == "running" ]]; then echo "failed" > "$CURRENT_STEP_FILE"; fi
  fi
}
on_signal() { echo "Interrupted"; mark_running_failed; exit 2; }
trap on_signal INT TERM HUP QUIT
trap mark_running_failed EXIT

# Defaults
coverage=10
threads=1
ignore_het="no"
run_craq="yes"
run_inspector="yes"
run_ragtag_corr="yes"
preset="standard"
resume=false
reads=""
craq_params=""
inspector_params=""
inspector_correct_params=""
ragtag_correct_params=""
ragtag_scf_params=""
# step5 (gap closing)
run_gap_closing="no"
tgsgapcloser_params=""

# Flag to force the raw preset for RagTag correct if no reads are given
force_raw_ragtag=false

# --- Helpers ----------------------------------------------------
q() { printf '%q' "$1"; }

# Absolute path resolver (no external deps required)
abs_path() {
  local p="$1"
  if [[ "$p" = /* ]]; then echo "$p"
  else
    local dir="$(cd "$(dirname "$p")" 2>/dev/null && pwd)"
    echo "$dir/$(basename "$p")"
  fi
}

# Build a final command string from a base array of tokens by:
#  1) removing any default options that appear in the user override string
#  2) appending the override string verbatim
# Each spec is "canon|syn1|syn2:has_value" (has_value=1 if the option expects a value).
_build_cmd_with_overrides() {
  local -n _base="$1"; local _override="$2"; shift 2
  local specs=("$@")

  declare -A to_remove=()
  for spec in "${specs[@]}"; do
    local keys="${spec%:*}"; local hasval="${spec##*:}"
    IFS='|' read -r -a names <<< "$keys"
    local present=false
    for nm in "${names[@]}"; do
      if [[ " ${_override} " =~ (^|[[:space:]])${nm}([=[:space:]]|$) ]]; then
        present=true; break
      fi
    done
    if $present; then
      local canon="${names[0]}"
      to_remove["$canon"]="$hasval"
      for nm in "${names[@]}"; do to_remove["$nm"]="$hasval"; done
    fi
  done

  local filtered=()
  local i=0
  while [[ $i -lt ${#_base[@]} ]]; do
    local tok="${_base[$i]}"; local skip=false; local consumes=0
    if [[ "${to_remove[$tok]+x}" == "x" ]]; then
      skip=true; consumes="${to_remove[$tok]}"
    fi
    if $skip; then
      if [[ "$consumes" == "1" && $((i+1)) -lt ${#_base[@]} ]]; then
        if [[ ! "${_base[$((i+1))]}" =~ ^- ]]; then ((i++)); fi
      fi
    else
      filtered+=("$tok")
    fi
    ((i++))
  done

  local out=""
  for t in "${filtered[@]}"; do out+=" $(q "$t")"; done
  out="${out# }"
  if [[ -n "$_override" ]]; then out+=" ${_override}"; fi
  printf '%s' "$out"
}

# --- Step-completion log checks------------------------
check_step_logs() {
  local step="$1"
  case "$step" in
    1)
      local f="$d1/craq.log"
      if [[ -f "$f" ]] && grep -Fqs "CRAQ analysis is finished." "$f" && ! grep -Fqs "Failed" "$f"; then
        echo "ok" > "$CURRENT_STEP_FILE"
      else
        echo "failed" > "$CURRENT_STEP_FILE"
      fi
      ;;
    2)
      local f1="$d2/Inspector_outputs/Inspector.log"
      local f2="$d2/Inspector_outputs/Inspector_correct.log"
      if [[ -f "$f1" && -f "$f2" ]] && grep -Fqs "Inspector evaluation finished. Bye." "$f1" && grep -Fqs "Inspector error correction finished. Bye." "$f2"; then
        echo "ok" > "$CURRENT_STEP_FILE"
      else
        echo "failed" > "$CURRENT_STEP_FILE"
      fi
      ;;
    3)
      local f="$d3/ragtag.correct.log"
      if [[ -f "$f" ]] && grep -Fqs "Goodbye" "$f"; then
        echo "ok" > "$CURRENT_STEP_FILE"
      else
        echo "failed" > "$CURRENT_STEP_FILE"
      fi
      ;;
    4)
      local f="$d4/ragtag.scaffold.log"
      if [[ -f "$f" ]] && grep -Fqs "Goodbye" "$f"; then
        echo "ok" > "$CURRENT_STEP_FILE"
      else
        echo "failed" > "$CURRENT_STEP_FILE"
      fi
      ;;
    5)
      local f="$d5/tgsgapcloser.log"
      if [[ -f "$f" ]] && grep -Fqs "ALL DONE !!!" "$f"; then
        echo "ok" > "$CURRENT_STEP_FILE"
      else
        echo "failed" > "$CURRENT_STEP_FILE"
      fi
      ;;
  esac
}
# -----------------------------------------------------------------------------#
#                                  Parse args                                  #
# -----------------------------------------------------------------------------#
if [[ $# -eq 0 ]]; then
  cat <<'EOF'
nohic-asm v1.0.0 - Contig correction and scaffolding pipeline
_____________________________________________________________

Usage: ./nohic-asm.sh -c <contigs.fa> -r <ref.fa> -o <outdir> [options]
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
EOF
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--contigs) contigs="$2"; shift 2 ;;
    -fq|--reads) reads="${2:-}"; shift 2 ;;
    -r|--reference) reference="$2"; shift 2 ;;
    -o|--output) outdir="$2"; shift 2 ;;
    -cov|--coverage) coverage="$2"; shift 2 ;;
    -t|--threads) threads="$2"; shift 2 ;;
    --ignore-het) ignore_het="$2"; shift 2 ;;
    --run-craq) run_craq="$2"; shift 2 ;;
    --run-inspector) run_inspector="$2"; shift 2 ;;
    --run-ragtag-correct) run_ragtag_corr="$2"; shift 2 ;;
    --run-gap-closing) run_gap_closing="$2"; shift 2 ;;
    -p|--presets) preset="$2"; shift 2 ;;
    --craq-params) craq_params="$2"; shift 2 ;;
    --inspector-params) inspector_params="$2"; shift 2 ;;
    --inspector-correct-params) inspector_correct_params="$2"; shift 2 ;;
    --ragtag-correct-params) ragtag_correct_params="$2"; shift 2 ;;
    --ragtag-scf-params) ragtag_scf_params="$2"; shift 2 ;;
    --tgsgapcloser-params) tgsgapcloser_params="$2"; shift 2 ;;
    --resume) resume=true; shift ;;
    -h|--help) "$0"; exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# Validation & normalization
if [[ -z "${contigs:-}" || -z "${reference:-}" || -z "${outdir:-}" ]]; then
  echo "Error: --contigs, --reference, and --output are required."; exit 1
fi
for v in run_craq run_inspector run_ragtag_corr run_gap_closing ignore_het; do
  val="$(echo "${!v}" | tr '[:upper:]' '[:lower:]')"
  if [[ "$v" == "ignore_het" ]]; then
    [[ "$val" =~ ^(yes|no)$ ]] || { echo "Error: --ignore-het must be yes or no"; exit 1; }
  else
    [[ "$val" =~ ^(yes|no)$ ]] || { echo "Error: --$v must be yes or no"; exit 1; }
  fi
  eval "$v=\"$val\""
done
[[ "$preset" =~ ^(standard|draft|luck|aggressive|raw)$ ]] || { echo "Error: --presets must be draft|luck|standard|aggressive|raw"; exit 1; }

# If no reads: disable CRAQ/Inspector; RagTag correct will run with 'raw' unless disabled
if [[ -z "$reads" ]]; then
  if [[ "$run_craq" == "yes" || "$run_inspector" == "yes" ]]; then
    echo "Warning: --reads not provided. Disabling CRAQ and Inspector; RagTag correct will run with 'raw' preset."
    run_craq="no"; run_inspector="no"
  fi
  if [[ "$run_ragtag_corr" != "no" ]]; then
    force_raw_ragtag=true
    run_ragtag_corr="yes"
  fi
fi

platform="unknown"
if [[ -n "$reads" ]]; then
  first_line=$( (gzip -dc "$reads" 2>/dev/null || head -n1 "$reads" 2>/dev/null || true) | head -n1 || true )
  if [[ "$first_line" == *@*runid=* ]]; then platform="ont"
  elif echo "$first_line" | grep -qi "ccs"; then platform="hifi"
  else platform="hifi"; fi
fi

if [[ "$platform" == "ont" ]]; then
  craq_map="map-ont"; ins_datatype_eval="nanopore"; ins_datatype_corr="nano-raw"; ragtag_read_type="ont"
else
  craq_map="map-hifi"; ins_datatype_eval="hifi"; ins_datatype_corr="pacbio-hifi"; ragtag_read_type="corr"
fi

# Layout & logging
mkdir -p "$outdir"
pipeline_log="$outdir/pipeline.log"
d1="$outdir/1_CRAQ"; d2="$outdir/2_Inspector"; d3="$outdir/3_RagTag_correct"; d4="$outdir/4_Scaffolding"; d5="$outdir/5_Gap_closing"
mkdir -p "$d1" "$d2" "$d3" "$d4" "$d5"

echo "nohic-asm started: $(date)" > "$pipeline_log"
echo "Threads=$threads Coverage=$coverage Platform=$platform" >> "$pipeline_log"
echo "Flags: CRAQ=$run_craq Inspector=$run_inspector RagTagCorrect=$run_ragtag_corr GapClosing=$run_gap_closing Preset=$preset Resume=$resume" >> "$pipeline_log"

step_status() { local f="$1"; [[ -f "$f" ]] || { echo "missing"; return; }; s="$(tr -d ' \t\r\n' < "$f" | tr '[:upper:]' '[:lower:]')"; [[ "$s" == "ok" ]] && echo "ok" || echo "failed"; }

start_step=1
if $resume; then
  s1=$(step_status "$d1/step_1_done.txt"); s2=$(step_status "$d2/step_2_done.txt"); s3=$(step_status "$d3/step_3_done.txt"); s4=$(step_status "$d4/step_4_done.txt"); s5=$(step_status "$d5/step_5_done.txt")

  # apply run-* flags when resuming: disabled steps are considered ok
  if [[ "$run_craq" == "no" ]]; then echo "ok" > "$d1/step_1_done.txt"; s1="ok"; fi
  if [[ "$run_inspector" == "no" ]]; then echo "ok" > "$d2/step_2_done.txt"; s2="ok"; fi
  if [[ "$run_ragtag_corr" == "no" ]]; then echo "ok" > "$d3/step_3_done.txt"; s3="ok"; fi
  if [[ "$run_gap_closing" == "no" ]]; then echo "ok" > "$d5/step_5_done.txt"; s5="ok"; fi

  if   [[ "$s1" != "ok" ]]; then start_step=1
  elif [[ "$s2" != "ok" ]]; then start_step=2
  elif [[ "$s3" != "ok" ]]; then start_step=3
  elif [[ "$s4" != "ok" ]]; then start_step=4
  elif [[ "$s5" != "ok" ]]; then start_step=5
  else echo "All steps already completed."; exit 0; fi

  echo "Resuming from step $start_step" | tee -a "$pipeline_log"
  for s in $(seq "$start_step" 5); do dir_var="d$s"; dir="${!dir_var}"; find "$dir" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true; done
else
  for dir in "$d1" "$d2" "$d3" "$d4" "$d5"; do find "$dir" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true; done
  start_step=1
fi

pick_craq_out()   { ls "$d1"/*.craq.* 2>/dev/null | head -n1 || ls "$d1"/*.craq 2>/dev/null | head -n1 || true; }
pick_insp_out()   { ls "$d2"/*.inspector.* 2>/dev/null | head -n1 || ls "$d2"/*.inspector 2>/dev/null | head -n1 || true; }
pick_rtcorr_out() { ls "$d3"/*.corrected.* 2>/dev/null | head -n1 || ls "$d3"/*.corrected 2>/dev/null | head -n1 || true; }
pick_scf_out()    { ls "$d4"/*.scf.* 2>/dev/null | head -n1 || ls "$d4"/*.scf 2>/dev/null | head -n1 || true; }

# --------------------------- Step 1: CRAQ (optional) ---------------------------
if [[ $start_step -le 1 ]]; then
  CURRENT_STEP_FILE="$d1/step_1_done.txt"
  if [[ "$run_craq" == "yes" ]]; then
    echo "running" > "$CURRENT_STEP_FILE"
    echo "[Step 1] Running CRAQ..." | tee -a "$pipeline_log"
    clip_rate="0.75"; [[ "$ignore_het" == "yes" ]] && clip_rate="0.55"

    craq_base=( craq --genome "$contigs" --sms_input "$reads" --break T --map "$craq_map" --output_dir "$d1/All_CRAQ_outputs"
                --sms_clip_coverRate "$clip_rate" --sms_coverage "$coverage" --thread "$threads" )
    
    craq_specs=( "--sms_coverage|-avgl:1" "--sms_clip_coverRate|-lf:1" "--thread|--threads|-t:1" "--map|-x:1" "--break:1" )

    craq_cmd=$(_build_cmd_with_overrides craq_base "$craq_params" "${craq_specs[@]}")
    echo "[CMD] $craq_cmd" >> "$pipeline_log"

    if eval "$craq_cmd" &> "$d1/craq.log"; then
      out="$d1/All_CRAQ_outputs/runAQI_out/out_correct.fa"
      if [[ -f "$out" ]]; then
        base=$(basename "$contigs"); b="${base%.*}"; e="${base##*.}"
        [[ "$b" == "$base" ]] && new="$d1/${base}.craq" || new="$d1/${b}.craq.${e}"
        mv "$out" "$new"
      fi
      echo "ok" > "$CURRENT_STEP_FILE"

      # verify logs -> set step status
      check_step_logs 1
      if [[ "$(tr -d ' \t\r\n' < "$CURRENT_STEP_FILE")" != "ok" ]]; then exit 1; fi    else
      echo "CRAQ failed (see $d1/craq.log)" | tee -a "$pipeline_log"; echo "failed" > "$CURRENT_STEP_FILE"; exit 1
    fi
  else
    echo "ok" > "$CURRENT_STEP_FILE"
  fi
fi

# ------------------------ Step 2: Inspector (optional) -------------------------
if [[ $start_step -le 2 ]]; then
  CURRENT_STEP_FILE="$d2/step_2_done.txt"
  if [[ "$run_inspector" == "yes" ]]; then
    echo "running" > "$CURRENT_STEP_FILE"
    echo "[Step 2] Running Inspector..." | tee -a "$pipeline_log"
    #mkdir -p "$d2/Inspector_outputs"
    insp_in="$(pick_craq_out)"; [[ -z "$insp_in" ]] && insp_in="$contigs"

    
    insp_eval_base=( inspector.py -c "$insp_in" -r "$reads" -o "$d2/Inspector_outputs"
                     -t "$threads" --datatype "$ins_datatype_eval" --min_contig_length 1000 )
    insp_eval_specs=( "-t|--thread|--threads:1" "--datatype|-d:1" "--min_contig_length:1" )

    insp_eval_cmd=$(_build_cmd_with_overrides insp_eval_base "$inspector_params" "${insp_eval_specs[@]}")
    echo "[CMD] $insp_eval_cmd" >> "$pipeline_log"

    if eval "$insp_eval_cmd" &> "$d2/inspector.log"; then
      
      insp_corr_base=( inspector-correct.py -i "$d2/Inspector_outputs" --datatype "$ins_datatype_corr" -t "$threads" -o "$d2/Inspector_outputs" )
      insp_corr_specs=( "-t|--thread|--threads:1" "--datatype:1" )

      insp_corr_cmd=$(_build_cmd_with_overrides insp_corr_base "$inspector_correct_params" "${insp_corr_specs[@]}")
      echo "[CMD] $insp_corr_cmd" >> "$pipeline_log"

      if eval "$insp_corr_cmd" &>> "$d2/inspector.log"; then
        corrected="$d2/Inspector_outputs/contig_corrected.fa"
        if [[ -f "$corrected" ]]; then
          in_name=$(basename "$insp_in"); b="${in_name%.*}"; e="${in_name##*.}"
          [[ "$b" == "$in_name" ]] && new="$d2/${in_name}.inspector" || new="$d2/${b}.inspector.${e}"
          mv "$corrected" "$new"
        fi
        echo "ok" > "$CURRENT_STEP_FILE"

      
      check_step_logs 2
      if [[ "$(tr -d ' \t\r\n' < "$CURRENT_STEP_FILE")" != "ok" ]]; then exit 1; fi      else
        echo "Inspector correction failed (see $d2/inspector.log)" | tee -a "$pipeline_log"; echo "failed" > "$CURRENT_STEP_FILE"; exit 1
      fi
    else
      echo "Inspector evaluation failed (see $d2/inspector.log)" | tee -a "$pipeline_log"; echo "failed" > "$CURRENT_STEP_FILE"; exit 1
    fi
  else
    echo "ok" > "$CURRENT_STEP_FILE"
  fi
fi

# ------------------- Step 3: RagTag correct (optional) -------------------------
if [[ $start_step -le 3 ]]; then
  CURRENT_STEP_FILE="$d3/step_3_done.txt"
  if [[ "$run_ragtag_corr" == "yes" ]]; then
    echo "running" > "$CURRENT_STEP_FILE"
    echo "[Step 3] Running RagTag correct..." | tee -a "$pipeline_log"

    rt_in="$(pick_insp_out)"; [[ -z "$rt_in" ]] && rt_in="$(pick_craq_out)"; [[ -z "$rt_in" ]] && rt_in="$contigs"

    # Build base from preset/defaults, then merge overrides
    rt_base=( ragtag.py correct "$reference" "$rt_in" -o "$d3" -t "$threads" )

    # Preset-driven defaults
    aligner="nucmer"; aligner_params=( --nucmer-params "--maxmatch -l 100 -c 500 -t $threads" )
    win=45000; remove_small="--remove-small"; extra=(); use_v=true; add_reads=true

    # If no reads were provided and RagTag correct is enabled, force 'raw'
    rt_preset="$preset"
    $force_raw_ragtag && rt_preset="raw"

    case "$rt_preset" in
      draft) aligner="minimap2"; aligner_params=( --mm2-params "-x asm5 -t $threads" ); win=10000; remove_small="";;
      luck)  aligner="minimap2"; aligner_params=( --mm2-params "-x asm5 -t $threads" ); win=45000; remove_small="--remove-small";;
      standard) : ;;
      aggressive) extra=( -d 50000 ); echo "Warning: aggressive RagTag correct selected." | tee -a "$pipeline_log" ;;
      raw)
        
        aligner="nucmer"; aligner_params=( --nucmer-params "--maxmatch -l 100 -c 500 -t $threads" )
        remove_small=""; use_v=false; extra=(); add_reads=false
        ;;
    esac
    # Add read validation only if allowed by preset and reads were provided
    if $add_reads && [[ -n "$reads" ]]; then
      rt_base+=( -R "$reads" -T "$ragtag_read_type" )
    fi

    rt_base+=( --aligner "$aligner" )
    rt_base+=( "${aligner_params[@]}" )
    $use_v && rt_base+=( -v "$win" )
    [[ -n "$remove_small" ]] && rt_base+=( "$remove_small" )
    [[ ${#extra[@]} -gt 0 ]] && rt_base+=( "${extra[@]}" )

    rt_specs=( "-t|--threads:1" "-R:1" "-T:1" "--aligner:1" "--mm2-params:1" "--nucmer-params:1" "-v:1" "-d:1" "--remove-small:0" )
    rt_cmd=$(_build_cmd_with_overrides rt_base "$ragtag_correct_params" "${rt_specs[@]}")

    echo "[CMD] $rt_cmd" >> "$pipeline_log"
    if eval "$rt_cmd" &> "$d3/ragtag.correct.log"; then
      out="$d3/ragtag.correct.fasta"
      if [[ -f "$out" ]]; then
        in_name=$(basename "$rt_in"); b="${in_name%.*}"; e="${in_name##*.}"
        
        if [[ "$b" == "$in_name" ]]; then
          new="$d3/${in_name}.corrected"
        else
          new="$d3/${b}.corrected.${e}"
        fi
        mv "$out" "$new"
      fi
      mkdir -p "$d3/RagTag_correct_outputs"
      
      for f in "$d3"/*; do
        fb=$(basename "$f")
        [[ "$fb" =~ ^(step_3_done\.txt|ragtag\.correct\.log|$(basename "${new:-x}")|RagTag_correct_outputs)$ ]] || mv "$f" "$d3/RagTag_correct_outputs/" 2>/dev/null || true
      done
      echo "ok" > "$CURRENT_STEP_FILE"

      
      check_step_logs 3
      if [[ "$(tr -d ' \t\r\n' < "$CURRENT_STEP_FILE")" != "ok" ]]; then exit 1; fi    else
      echo "RagTag correct failed (see $d3/ragtag.correct.log)" | tee -a "$pipeline_log"; echo "failed" > "$CURRENT_STEP_FILE"; exit 1
    fi
  else
    echo "ok" > "$CURRENT_STEP_FILE"
  fi
fi

# --------------------- Step 4: RagTag scaffold (required) ----------------------
if [[ $start_step -le 4 ]]; then
  CURRENT_STEP_FILE="$d4/step_4_done.txt"
  echo "running" > "$CURRENT_STEP_FILE"
  echo "[Step 4] Running RagTag scaffold..." | tee -a "$pipeline_log"

  sc_in="$(pick_rtcorr_out)"; [[ -z "$sc_in" ]] && sc_in="$(pick_insp_out)"; [[ -z "$sc_in" ]] && sc_in="$(pick_craq_out)"; [[ -z "$sc_in" ]] && sc_in="$contigs"
  sc_cmd="ragtag.py scaffold $(q "$reference") $(q "$sc_in") -o $(q "$d4") -t $(q "$threads")"
  # append user-specified RagTag scaffold params if provided
  if [[ -n "$ragtag_scf_params" ]]; then sc_cmd+=" $ragtag_scf_params"; fi
  echo "[CMD] $sc_cmd" >> "$pipeline_log"

  if eval "$sc_cmd" &> "$d4/ragtag.scaffold.log"; then
    out="$d4/ragtag.scaffold.fasta"
    if [[ -f "$out" ]]; then
      in_name=$(basename "$sc_in"); b="${in_name%.*}"; e="${in_name##*.}"
      [[ "$b" == "$in_name" ]] && new="$d4/${in_name}.scf" || new="$d4/${b}.scf.${e}"
      mv "$out" "$new"
    fi
    mkdir -p "$d4/RagTag_scaffold_outputs"
    for f in "$d4"/*; do
      fb=$(basename "$f")
      [[ "$fb" =~ ^(step_4_done\.txt|ragtag\.scaffold\.log|$(basename "${new:-x}")|RagTag_scaffold_outputs)$ ]] || mv "$f" "$d4/RagTag_scaffold_outputs/" 2>/dev/null || true
    done
    echo "ok" > "$CURRENT_STEP_FILE"

      
      check_step_logs 4
      if [[ "$(tr -d ' \t\r\n' < "$CURRENT_STEP_FILE")" != "ok" ]]; then exit 1; fi  else
    echo "RagTag scaffold failed (see $d4/ragtag.scaffold.log)" | tee -a "$pipeline_log"; echo "failed" > "$CURRENT_STEP_FILE"; exit 1
  fi
fi

# --------------------- Step 5: TGSGapCloser (optional) -------------------------
if [[ $start_step -le 5 ]]; then
  CURRENT_STEP_FILE="$d5/step_5_done.txt"
  if [[ "$run_gap_closing" == "yes" ]]; then
    echo "running" > "$CURRENT_STEP_FILE"
    echo "[Step 5] Running TGSGapCloser..." | tee -a "$pipeline_log"

    scf_src="$(pick_scf_out)"
    if [[ -z "$scf_src" || ! -f "$scf_src" ]]; then
      echo "Error: no scaffolded assembly (*.scf.*) found in $d4" | tee -a "$pipeline_log"
      echo "failed" > "$CURRENT_STEP_FILE"; exit 1
    fi
    if [[ -z "$reads" ]]; then
      echo "Error: --reads is required for TGSGapCloser (needs FASTA reads). Provide -fq/--reads." | tee -a "$pipeline_log"
      echo "failed" > "$CURRENT_STEP_FILE"; exit 1
    fi

    # Resolve absolute paths
    abs_scf="$(abs_path "$scf_src")"
    abs_reads="$(abs_path "$reads")"

    # Convert FASTQ->FASTA
    ( cd "$d5" && seqkit fq2fa -j "$threads" "$abs_reads" -o reads.fa ) >> "$pipeline_log" 2>&1

    # Build TGSGapCloser command
    scf_basename="$(basename "$abs_scf")"
    scf_noext="${scf_basename%.*}"
    tgs_prefix="${scf_noext}.tgs"
    tgs_type="$([[ "$platform" == "ont" ]] && echo "ont" || echo "pb")"

    tgs_base=( tgsgapcloser
               --scaff "$abs_scf"
               --reads "reads.fa"
               --output "$tgs_prefix"
               --minmap_arg "-x asm5 -t $threads"
               --ne
               --tgstype "$tgs_type"
               --thread "$threads" )

    # Overridable specs
    tgs_specs=( "--scaff:1" "--reads:1" "--output:1" "--minmap_arg:1" "--tgstype:1" "--thread|-t:1" )

    tgs_cmd=$(_build_cmd_with_overrides tgs_base "$tgsgapcloser_params" "${tgs_specs[@]}")

    echo "[CMD] (cd $d5 && $tgs_cmd)" >> "$pipeline_log"
    if ( cd "$d5" && eval "$tgs_cmd" ) &> "$d5/tgsgapcloser.log"; then
      # Rename <prefix>.scaff_seq -> <prefix>.fa
      final_src="$d5/${tgs_prefix}.scaff_seqs"
      final_dst="$d5/${tgs_prefix}.fa"
      if [[ -f "$final_src" ]]; then
        mv "$final_src" "$final_dst"
      fi
      mkdir -p "$d5/TGSGapcloser_outputs"
      
      for f in "$d5"/*; do
        fb=$(basename "$f")
        if [[ "$fb" == "$(basename "$final_dst")" || "$fb" == "TGSGapcloser_outputs" || "$fb" == "step_5_done.txt" || "$fb" == "tgsgapcloser.log" ]]; then
          continue
        fi
        mv "$f" "$d5/TGSGapcloser_outputs/" 2>/dev/null || true
      done
      echo "ok" > "$CURRENT_STEP_FILE"

      
      check_step_logs 5
      if [[ "$(tr -d ' \t\r\n' < "$CURRENT_STEP_FILE")" != "ok" ]]; then exit 1; fi    else
      echo "TGSGapCloser failed (see $d5/tgsgapcloser.log)" | tee -a "$pipeline_log"; echo "failed" > "$CURRENT_STEP_FILE"; exit 1
    fi
  else
    echo "ok" > "$CURRENT_STEP_FILE"
  fi
fi

echo "nohic-asm was done at: $(date)" | tee -a "$pipeline_log"

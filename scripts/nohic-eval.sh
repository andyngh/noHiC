#!/usr/bin/env bash
# nohic-eval.sh v1.0.0 — Assembly evaluation pipeline
set -euo pipefail

# --------------------------- helpers ---------------------------

version="1.0.0"

log() {
  local ts
  ts="$(date '+%Y-%m-%d %H:%M:%S')"
  echo "[$ts] $*" | tee -a "${PIPELINE_LOG}"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

# marker writer: $1=dir, $2=filename, $3=ok|failed
write_marker() {
  local dir="$1"; local fname="$2"; local val="$3"
  mkdir -p "$dir"
  echo "$val" > "${dir}/${fname}"
}

# check non-empty file(s)
nonempty() {
  for f in "$@"; do
    [[ -s "$f" ]] || return 1
  done
  return 0
}

# normalize absolute path
abspath() {
  python3 - "$1" <<'PY'
import os, sys; print(os.path.abspath(sys.argv[1]))
PY
}

current_step_dir=""
current_step_marker=""
on_interrupt() {
  if [[ -n "${current_step_dir}" && -n "${current_step_marker}" ]]; then
    write_marker "${current_step_dir}" "${current_step_marker}" "failed"
  fi
  log "Received termination signal. Marked ${current_step_marker} as failed (if set). Exiting."
  exit 130
}
on_error() {
  local ec=$?
  if [[ -n "${current_step_dir}" && -n "${current_step_marker}" ]]; then
    write_marker "${current_step_dir}" "${current_step_marker}" "failed"
  fi
  log "Step failed with exit code ${ec}. Marked ${current_step_marker} as failed (if set)."
  exit "${ec}"
}
trap on_interrupt INT TERM HUP
trap on_error ERR

# --------------------------- defaults ---------------------------

RUN_BUSCO="yes"
RUN_CRAQ="yes"
RUN_INSPECTOR="yes"
RUN_VIZ="yes"
RESUME="no"
BUSCO_DB=""
THREADS=""
OUTDIR=""
ASSEMBLY=""
READS=""
COVERAGE=""
PLATFORM=""
CHR_NAME_FILE=""
MM2_PARAMS="-cx asm10"
REFERENCE=""
SCALE_FACTOR="500000"

# --------------------------- usage ---------------------------

print_help() {
cat <<'EOF'
nohic-eval.sh v1.0.0 — Assembly evaluation pipeline
_________________________________________________

Required:
      --input-assembly, -i <.fasta|.fa|.fna>    Assembly file (can be gzipped)
      --output-directory, -o <dir>    		      Common output directory containing all evaluation results
      --threads, -t <num>             		      Number of threads

Optional evaluations (default: yes):
      --run-busco <yes|no>        		          Run BUSCO
      --run-craq <yes|no>         		          Run CRAQ
      --run-inspector <yes|no>    		          Run Inspector
      --visualization, -v <yes|no>    		      Visualize assembly errors

Other optional args:
      --busco-db, -b <str>           		        BUSCO lineage database to be used (required if --run-busco yes)
      --reads, -r <fasq[.gz]>              	    Long reads for CRAQ/Inspector
      --coverage <num>            		          Sequencing coverage for CRAQ
      --sequencing-platform, -p <pb|hifi|ont>	  Sequencing platform
      --chr-name <.txt>           		          Text file of chromosome names (one per line) for visualization
      --scale-factor <int>			                Stretch misassemblies to this length (bp) for visualization [default: 500000]
      --mm2-params "<str>"               	      Minimap2 params for dot plot [default: "-cx asm10"]
      --reference <fasta>                 	    Reference genome FASTA for dot plot
      --resume                  		            Resume nohic-eval from the earliest failed/missing step
      --help, -h                    		        Show this help
      --version                     		        Show version

EOF
}

# --------------------------- arg parsing ---------------------------

long_opts="input-assembly:,output-directory:,threads:,run-busco:,run-craq:,run-inspector:,visualization:,busco-db:,reads:,coverage:,sequencing-platform:,chr-name:,mm2-params:,reference:,resume,help,version"
short_opts="i:o:t:b:r:p:v:h"

# Use getopt if available for long options
if getopt --test >/dev/null 2>&1; then
  PARSED=$(getopt --options="$short_opts" --longoptions="$long_opts" --name "nohic-eval" -- "$@") || { print_help; exit 2; }
  eval set -- "$PARSED"
fi

while true; do
  case "${1-}" in
    -i|--input-assembly) ASSEMBLY="$2"; shift 2;;
    -o|--output-directory) OUTDIR="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    --run-busco) RUN_BUSCO="$2"; shift 2;;
    --run-craq) RUN_CRAQ="$2"; shift 2;;
    --run-inspector) RUN_INSPECTOR="$2"; shift 2;;
    -v|--visualization) RUN_VIZ="$2"; shift 2;;
    -b|--busco-db) BUSCO_DB="$2"; shift 2;;
    -r|--reads) READS="$2"; shift 2;;
    --coverage) COVERAGE="$2"; shift 2;;
    -p|--sequencing-platform) PLATFORM="$2"; shift 2;;
    --chr-name) CHR_NAME_FILE="$2"; shift 2;;
    --scale-factor) SCALE_FACTOR="$2"; shift 2;;
    --mm2-params) MM2_PARAMS="$2"; shift 2;;
    --reference) REFERENCE="$2"; shift 2;;
    --resume) RESUME="yes"; shift ;;
    -h|--help) print_help; exit 0;;
    --version) echo "$version"; exit 0;;
    --) shift; break;;
    *) break;;
  esac
done

# validate mandatory
[[ -n "$ASSEMBLY" ]] || { print_help; die "Missing --input-assembly"; }
[[ -n "$OUTDIR"   ]] || { print_help; die "Missing --output-directory"; }
[[ -n "$THREADS"  ]] || { print_help; die "Missing --threads"; }

# resolve OUTDIR to absolute to avoid path issues after `cd`
mkdir -p "$OUTDIR"
OUTDIR_ABS="$(abspath "$OUTDIR")"

PIPELINE_LOG="${OUTDIR_ABS}/pipeline.log"
: > "$PIPELINE_LOG"

# record pipeline start time for total duration report
PIPELINE_START_TS="$(date +%s)"

log "nohic-eval.sh v${version}"
ASSEMBLY_ABS="$(abspath "$ASSEMBLY")"
ASSEMBLY_NAME="$(basename "$ASSEMBLY_ABS")"
log "Input assembly: $ASSEMBLY_ABS"
log "Outdir: $OUTDIR_ABS"
log "Threads: $THREADS"
log "Flags: busco=$RUN_BUSCO craq=$RUN_CRAQ inspector=$RUN_INSPECTOR viz=$RUN_VIZ resume=$RESUME"
log "Reads: ${READS:-<none>} | Coverage: ${COVERAGE:-<none>} | Platform: ${PLATFORM:-<none>} | BUSCO DB: ${BUSCO_DB:-<none>}"
[[ -n "$CHR_NAME_FILE" ]] && log "Chromosome name file: $CHR_NAME_FILE"

log "Reference genome for drawing dot plot: ${REFERENCE:-<none>}"
log "minimap2 params for drawing dot plot: ${MM2_PARAMS}"

# create subdirs (absolute)
D1="${OUTDIR_ABS}/1_Assembly_statistics"
D2="${OUTDIR_ABS}/2_BUSCO"
D3="${OUTDIR_ABS}/3_CRAQ"
D4="${OUTDIR_ABS}/4_Inspector"
D5="${OUTDIR_ABS}/5_Error_visualization"
mkdir -p "$D1" "$D2" "$D3" "$D4" "$D5"

# prepare working assembly (handle .gz) — keep original filename for all downstream tools
WORK="${OUTDIR_ABS}/tmp"
mkdir -p "$WORK"
WORK_ASM="${WORK}/${ASSEMBLY_NAME}"
if [[ "$ASSEMBLY_ABS" =~ \.gz$ ]]; then
  log "Decompressing gzipped assembly to $WORK_ASM"
  gzip -cd "$ASSEMBLY_ABS" > "$WORK_ASM"
else
  log "Linking assembly to $WORK_ASM"
  ln -sf "$ASSEMBLY_ABS" "$WORK_ASM"
fi

# derive mapping preset / datatype from platform
MAP_PRESET=""
SEQ_DATATYPE=""
case "${PLATFORM:-}" in
  pb)   MAP_PRESET="map-pb";   SEQ_DATATYPE="clr" ;;
  hifi) MAP_PRESET="map-hifi"; SEQ_DATATYPE="hifi";;
  ont)  MAP_PRESET="map-ont";  SEQ_DATATYPE="nanopore";;
  "" )  : ;;
  * )   die "Invalid --sequencing-platform (use pb|hifi|ont)";;
esac

# --------------------------- step selection (resume) ---------------------------
# Determine starting step based on markers and toggles.
# step indices: 1 metrics, 2 busco, 3 craq, 4 inspector, 5 viz

read_marker() {
  local f="$1"
  [[ -s "$f" ]] && awk '{print $1}' "$f" || echo "missing"
}

need_step1="yes"
need_step2="$RUN_BUSCO"
need_step3="$RUN_CRAQ"
need_step4="$RUN_INSPECTOR"
need_step5="$RUN_VIZ"

# auto turn off step5 if both step3 and step4 are "no"
if [[ "$RUN_CRAQ" == "no" && "$RUN_INSPECTOR" == "no" ]]; then
  need_step5="no"
fi

start_step=1
if [[ "$RESUME" == "yes" ]]; then
  m1=$(read_marker "${D1}/step_1_done.txt")
  m2=$(read_marker "${D2}/step_2_done.txt")
  m3=$(read_marker "${D3}/step_3_done.txt")
  m4=$(read_marker "${D4}/step_4_done.txt")
  m5=$(read_marker "${D5}/step_5_done.txt")

  # helper to set earliest needed step
  start_step=6
  if [[ "$need_step1" == "yes" && "$m1" != "ok" ]]; then start_step=1; fi
  if [[ "$need_step2" == "yes" && "$m2" != "ok" && $start_step -gt 2 ]]; then start_step=2; fi
  if [[ "$need_step3" == "yes" && "$m3" != "ok" && $start_step -gt 3 ]]; then start_step=3; fi
  if [[ "$need_step4" == "yes" && "$m4" != "ok" && $start_step -gt 4 ]]; then start_step=4; fi
  if [[ "$need_step5" == "yes" && "$m5" != "ok" && $start_step -gt 5 ]]; then start_step=5; fi
  if [[ $start_step -eq 6 ]]; then
    log "All requested steps are already ok. Nothing to do."
    exit 0
  fi
  log "Resume: restarting from step ${start_step}"
else
  start_step=1
fi

# --------------------------- step 1: assembly metrics ---------------------------
if [[ $start_step -le 1 && "$need_step1" == "yes" ]]; then
  current_step_dir="$D1"; current_step_marker="step_1_done.txt"
  log "Step 1: Assembly metrics"
  # log main command (single line)
  log "COMMAND: gfastats \"$WORK_ASM\" -j \"$THREADS\" -t > \"${D1}/assembly_stats.txt\""
  (
    set -x
    gfastats "$WORK_ASM" -j "$THREADS" -t > "${D1}/assembly_stats.txt"
    bioawk -c fastx '{print $name, length($seq)}' < "$WORK_ASM" | sort -n -r -k2 > "${D1}/scaffold_lengths.txt"
  ) &>> "${D1}/stats.log"

  if nonempty "${D1}/assembly_stats.txt" "${D1}/scaffold_lengths.txt"; then
    write_marker "$D1" "$current_step_marker" "ok"
  else
    write_marker "$D1" "$current_step_marker" "failed"
  fi
fi

# --------------------------- step 2: BUSCO ---------------------------
if [[ $start_step -le 2 && "$need_step2" == "yes" ]]; then
  current_step_dir="$D2"; current_step_marker="step_2_done.txt"
  mkdir -p "$D2"
  log "Step 2: BUSCO (db=$BUSCO_DB)"
  : > "${D2}/busco.log"
  if [[ -z "$BUSCO_DB" ]]; then
    log "BUSCO requested but --busco-db not provided"
    write_marker "$D2" "$current_step_marker" "failed"
  else
    (
      set -x
      cd "$D2"
      busco --download "$BUSCO_DB"
      LINEAGE_DIR="${D2}/busco_downloads/lineages/${BUSCO_DB}"
      BUSCO_OUTNAME="${ASSEMBLY_NAME}.busco"
      # log main command (single line)
      log "COMMAND: busco -i \"$WORK_ASM\" -m genome --cpu \"$THREADS\" -l \"$LINEAGE_DIR\" --out_path \"$D2\" -o \"$BUSCO_OUTNAME\""
      busco -i "$WORK_ASM" -m genome --cpu "$THREADS" -l "$LINEAGE_DIR" \
            --out_path "$D2" -o "$BUSCO_OUTNAME" --skip_bbtools #Skipping BBTools to avoid OOM error
    ) &>> "${D2}/busco.log"

    # normalize run_<name> -> <name>
    if [[ -d "${D2}/run_${ASSEMBLY_NAME}.busco" && ! -e "${D2}/${ASSEMBLY_NAME}.busco" ]]; then
      mv -f "${D2}/run_${ASSEMBLY_NAME}.busco" "${D2}/${ASSEMBLY_NAME}.busco"
    fi

    BUSCO_OUT="${D2}/${ASSEMBLY_NAME}.busco"
    n_ss=$(find "$BUSCO_OUT" -maxdepth 1 -type f -name 'short_summary.*' | wc -l | awk '{print $1}')
    if [[ "$n_ss" -ge 2 ]]; then
      write_marker "$D2" "$current_step_marker" "ok"
    else
      write_marker "$D2" "$current_step_marker" "failed"
    fi
  fi
elif [[ "$need_step2" == "no" ]]; then
  write_marker "$D2" "step_2_done.txt" "ok"
fi

# --------------------------- step 3: CRAQ ---------------------------
if [[ $start_step -le 3 && "$need_step3" == "yes" ]]; then
  current_step_dir="$D3"; current_step_marker="step_3_done.txt"
  log "Step 3: CRAQ"
  : > "${D3}/craq.log"
  if [[ -z "$READS" || -z "$PLATFORM" ]]; then
    log "CRAQ requested but --reads and/or --sequencing-platform missing"
    write_marker "$D3" "$current_step_marker" "failed"
  else
    # log main command (single line)
    log "COMMAND: craq --genome \"$WORK_ASM\" --sms_input \"$READS\" --sms_coverage \"${COVERAGE:-30}\" --break F --map \"$MAP_PRESET\" --thread \"$THREADS\" --output_dir \"${D3}/All_CRAQ_outputs\""
    (
      set -x
      craq --genome "$WORK_ASM" \
           --sms_input "$READS" \
           --sms_coverage "${COVERAGE:-30}" \
           --break F \
           --map "$MAP_PRESET" \
           --thread "$THREADS" \
           --output_dir "${D3}/All_CRAQ_outputs"
    ) &>> "${D3}/craq.log"
    if grep -q "CRAQ analysis is finished." "${D3}/craq.log" && ! grep -q "Failed" "${D3}/craq.log"; then
      write_marker "$D3" "$current_step_marker" "ok"
    else
      write_marker "$D3" "$current_step_marker" "failed"
    fi
  fi
elif [[ "$need_step3" == "no" ]]; then
  write_marker "$D3" "step_3_done.txt" "ok"
fi

# post-processing for CRAQ (only if ok or step skipped==no? Spec: do tasks if step_3_done.txt has ok)
if [[ "$(read_marker "${D3}/step_3_done.txt")" == "ok" && "$need_step3" == "yes" ]]; then
  if [[ -f "${D3}/All_CRAQ_outputs/runAQI_out/out_final.Report" ]]; then
    mv -f "${D3}/All_CRAQ_outputs/runAQI_out/out_final.Report" "${D3}/CRAQ_AQI_metrics.txt" || true
  fi
  CSE_IN="${D3}/All_CRAQ_outputs/LRout/LR_putative.SE"
  CSE_OUT="${D3}/CSE.csv"
  if [[ -s "$CSE_IN" ]]; then
    awk 'BEGIN{OFS=","} {print $1,$2,$2+1,"CSE"}' "$CSE_IN" > "$CSE_OUT"
  else
    : > "$CSE_OUT"
  fi
else
    : > "${D3}/CSE.csv"
fi

# --------------------------- step 4: Inspector ---------------------------
if [[ $start_step -le 4 && "$need_step4" == "yes" ]]; then
  current_step_dir="$D4"; current_step_marker="step_4_done.txt"
  mkdir -p "${D4}/All_inspector_outputs"
  log "Step 4: Inspector"
  : > "${D4}/inspector.log"
  if [[ -z "$READS" || -z "$PLATFORM" ]]; then
    log "Inspector requested but --reads and/or --sequencing-platform missing"
    write_marker "$D4" "$current_step_marker" "failed"
  else
    # log main command (single line)
    log "COMMAND: inspector.py -c \"$WORK_ASM\" -r \"$READS\" -o \"${D4}/All_inspector_outputs\" -t \"$THREADS\" --datatype \"$SEQ_DATATYPE\""
    (
      set -x
      inspector.py -c "$WORK_ASM" \
                   -r "$READS" \
                   -o "${D4}/All_inspector_outputs" \
                   -t "$THREADS" \
                   --datatype "$SEQ_DATATYPE"
    ) &>> "${D4}/inspector.log"
    if grep -q "Inspector evaluation finished. Bye." "${D4}/All_inspector_outputs/Inspector.log"; then
      write_marker "$D4" "$current_step_marker" "ok"
    else
      write_marker "$D4" "$current_step_marker" "failed"
    fi
  fi
elif [[ "$need_step4" == "no" ]]; then
  write_marker "$D4" "step_4_done.txt" "ok"
fi

# post-processing for Inspector if ok
# (UPDATED)

if [[ "$(read_marker "${D4}/step_4_done.txt")" == "ok" && "$need_step4" == "yes" ]]; then
  # moving summary if present
  if [[ -f "${D4}/All_inspector_outputs/summary_statistics" ]]; then
    mv -f "${D4}/All_inspector_outputs/summary_statistics" "${D4}/Assembly_statistics_Inspector.txt" || true
  fi

  # 1) small_scale_error.bed -> small_scale_error.csv (strict: must exist)
  SML_IN="${D4}/All_inspector_outputs/small_scale_error.bed"
  SML_OUT="${D4}/small_scale_error.csv"
  [[ -f "$SML_IN" ]] || die "Missing ${SML_IN}"
  awk 'BEGIN{FS="\t"; OFS=","} NR>1 && $9<0.05 {print $1,$2,$3,$8}' "$SML_IN" > "$SML_OUT"

  # 2) structural_error.bed -> structural_error.csv (strict: must exist, drop header)
  STR_IN="${D4}/All_inspector_outputs/structural_error.bed"
  STR_OUT="${D4}/structural_error.csv"
  [[ -f "$STR_IN" ]] || die "Missing ${STR_IN}"
  awk 'BEGIN{FS="\t"; OFS=","} NR>1 {print $1,$2,$3,$5}' "$STR_IN" > "$STR_OUT"

  # 3) Merge into inspector_errors.csv with header
  INS_COMB="${D4}/inspector_errors.csv"
  {
    echo "chrom,start,end,type"
    cat "$STR_OUT" "$SML_OUT"
  } > "$INS_COMB"
else
    echo "chrom,start,end,type" > "${D4}/inspector_errors.csv"
fi


# --------------------------- step 5: Visualization ---------------------------
if [[ $start_step -le 5 && "$need_step5" == "yes" ]]; then
  current_step_dir="$D5"; current_step_marker="step_5_done.txt"
  log "Step 5: Visualization"
  ERR_COMB="${D5}/error_to_plots.csv"
  { 
    if [[ -s "${D4}/inspector_errors.csv" ]]; then
      cat "${D4}/inspector_errors.csv"
    else
      echo "chrom,start,end,type"
    fi
    if [[ -s "${D3}/CSE.csv" ]]; then
      cat "${D3}/CSE.csv"
    fi
  } > "$ERR_COMB"

  SCAFF="${D1}/scaffold_lengths.txt"
  CHR_LEN="${D5}/chr_len.csv"
  if [[ -s "$SCAFF" ]]; then
    if [[ -n "$CHR_NAME_FILE" && -s "$CHR_NAME_FILE" ]]; then
      echo "chrom,length" > "$CHR_LEN"
      awk '{print $1","$2}' "$SCAFF" > "${D5}/__all_chr_len.csv"
      while IFS= read -r nm; do
        grep -E "^${nm}," "${D5}/__all_chr_len.csv" || true
      done < "$CHR_NAME_FILE" >> "$CHR_LEN"
      rm -f "${D5}/__all_chr_len.csv"
    else
      { echo "chrom,length"; awk '{print $1","$2}' "$SCAFF"; } > "$CHR_LEN"
    fi
  else
    die "Missing scaffold_lengths.txt for visualization"
  fi

  VIZ_R="${D5}/nohic-viz.R"
  cat > "$VIZ_R" <<'RSCRIPT'
#!/usr/bin/env Rscript

ensure_pkgs <- function(pkgs) {
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (user_lib == "") {
    user_lib <- file.path(path.expand("~"), "R", "library")
    Sys.setenv(R_LIBS_USER = user_lib)
  }
  if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install)) {
    install.packages(to_install, repos = "https://cloud.r-project.org", lib = user_lib, dependencies = TRUE)
  }
}
ensure_pkgs(c("ggplot2", "readr", "dplyr", "svglite"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: nohic-viz.R <chrom_lengths.csv> <misassemblies.csv> <output_prefix> <scale_factor>")
}
chrom_csv <- args[1]
mis_csv   <- args[2]
out_pref  <- args[3]
scale_factor <- args[4]

chroms <- readr::read_csv(chrom_csv, show_col_types = FALSE, col_names = TRUE)
if (!all(c("chrom", "length") %in% names(chroms))) {
  if (ncol(chroms) >= 2) {
    names(chroms)[1:2] <- c("chrom","length")
  } else stop("First CSV must have two columns: chrom,length (with header).")
}
chroms <- chroms %>% mutate(chrom = trimws(as.character(chrom)))

mises <- readr::read_csv(mis_csv, show_col_types = FALSE, col_names = TRUE)
if (!all(c("chrom", "start", "end", "type") %in% names(mises))) {
  if (ncol(mises) >= 4) {
    names(mises)[1:4] <- c("chrom","start","end","type")
  } else stop("Second CSV must have three columns: chrom,start,end,type (with header).")
}
mises  <- mises  %>% mutate(chrom = trimws(as.character(chrom)))

mises <- mises %>% dplyr::semi_join(chroms, by = "chrom")
mises <- mises %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  mutate(x0 = pmin(start, end), x1 = pmax(start, end)) %>%
  select(chrom, x0, x1, type)
mises <- mises %>% mutate(type = as.factor(type))

chrom_lengths <- chroms %>% select(chrom, length)
mises <- mises %>%
  inner_join(chrom_lengths, by = "chrom") %>%
  mutate(x0 = pmax(0, pmin(x0, length)),
         x1 = pmax(0, pmin(x1, length))) %>%
  filter(x1 > x0)

#min_vis_bp <- max(5e4, 5e-5 * max(chroms$length, na.rm = TRUE))
min_vis_bp <- as.numeric(scale_factor)

mises <- mises %>%
  mutate(len = x1 - x0,
         d   = pmax(len, min_vis_bp),
         mid = (x0 + x1) / 2,
         x0_vis = pmax(0, pmin(length, mid - d/2)),
         x1_vis = pmax(0, pmin(length, mid + d/2)))

chroms$chrom <- factor(chroms$chrom, levels = chroms$chrom)
mises$chrom  <- factor(mises$chrom,  levels = levels(chroms$chrom))

chrom_seg <- chroms %>% transmute(chrom, x0 = 0, x1 = length)
line_mm <- 4

p <- ggplot() +
  geom_segment(data = chrom_seg,
               aes(x = x0, xend = x1, y = chrom, yend = chrom),
               linewidth = line_mm, color = "black", lineend = "round") +
  geom_segment(data = mises,
               aes(x = x0_vis, xend = x1_vis, y = chrom, yend = chrom, colour = type),
               linewidth = line_mm, lineend = "butt") +
  scale_colour_discrete(name = "Misassembly Types")+
  scale_y_discrete(name = NULL, drop = FALSE) +
  scale_x_continuous(name = "Position", expand = expansion(mult = c(0.02, 0.01))) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.ticks.y     = element_blank())

max_len <- max(chroms$length, na.rm = TRUE)
n_chr   <- nrow(chroms)
w_in <- max(6, min(16, max_len / 50e6 * 10))
h_in <- max(2.5, min(20, n_chr * 0.35))

ggsave(filename = paste0(out_pref, ".svg"), plot = p, width = w_in, height = h_in, units = "in", device = "svg")
ggsave(filename = paste0(out_pref, ".png"), plot = p, width = w_in, height = h_in, units = "in", dpi = 1200)
RSCRIPT
  chmod +x "$VIZ_R"

  # log main command (single line)
  log "COMMAND: \"$VIZ_R\" \"$CHR_LEN\" \"$ERR_COMB\" \"${D5}/plotted_errors\" \"$SCALE_FACTOR\""
  ( set -x; "$VIZ_R" "$CHR_LEN" "$ERR_COMB" "${D5}/plotted_errors" "$SCALE_FACTOR" ) &>> "${PIPELINE_LOG}"

  # ---------------- create paf2dotplot.R and run dot-plot ----------------
  PAF2_R="${D5}/paf2dotplot.R"
  cat > "$PAF2_R" <<'RS2'
#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))

generate_colors <- function(n){
  # check RColorBrewer
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required for generate_colors. Please install it with install.packages('RColorBrewer')")
  }

  if (n <= 12) {
    RColorBrewer::brewer.pal(n, "Set3")[1:n]
  } else {
    # When n > 12, extend using colorRampPalette
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
  }
}

option_list <- list(
  make_option(c("-o","--output"), type="character",
              help="output filename prefix [input.paf]", 
              dest="output_filename"),
  make_option(c("-p","--plot-size"), type="numeric", default=15,
              help="plot size X inches [%default]",
              dest="plot_size"),
  make_option(c("-f","--flip"), action="store_true", default=FALSE,
              help="flip query if most alignments are in reverse complement [%default]",
              dest="flip"),
  make_option(c("-b","--break-point"), action="store_true", default=FALSE,
              help="show break points [%default]",
              dest="break_point"),
  make_option(c("-c","--identity-lower-color"), type="numeric", default=0,
              help="percent identity that lower than this value will be assigned the same color [%default]",
              dest="identity_lower"),
  make_option(c("-s", "--sort-by-refid"), action="store_true", default=FALSE,
              help="sort reference IDs in alphabetical order, default by length [%default]",
              dest="sortbyID"),
  make_option(c("-q", "--min-query-length"), type="numeric", default=400000,
              help="filter queries with total alignments less than cutoff X bp [%default]",
              dest="min_query_aln"),
  make_option(c("-m", "--min-alignment-length"), type="numeric", default=10000,
              help="filter alignments less than cutoff X bp [%default]",
              dest="min_align"),
  make_option(c("-r", "--min-ref-len"), type="numeric", default=1000000,
              help="filter references with length less than cutoff X bp [%default]",
              dest="min_ref_len"),
  make_option(c("-i", "--reference-ids"), type="character", default=NULL,
              help="comma-separated list of reference IDs to keep and order [%default]",
              dest="refIDs"),
  make_option(c("-e", "--ref-bed"), type="character", default=NULL,
              help="reference BED file to draw vertical dotted lines at the specified regions.",
              dest="ref_bed_file"),
  make_option(c("-E", "--query-bed"), type="character", default=NULL,
              help="query BED file to draw horizontal dotted lines at the specified regions.",
              dest="query_bed_file"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Show version and exit")
)

options(error=traceback)
parser <- OptionParser(usage = "%prog [options] input.paf\n\nFor more information, see https://github.com/moold/paf2dotplot", option_list = option_list)
opts = parse_args(parser, positional_arguments = c(0, 1))
opt = opts$options

script_version <- "1.0.1"
if (opt$version) {
  cat(paste0("version: ", script_version, "\n"))
  quit(status=0)
}

input_file = opts$args
if(length(input_file) <= 0){
  cat(sprintf("Error: missing input file: input.paf!\n\n"))
  print_help(parser)
  quit()
}else if (file.access(input_file, mode=4) == -1){
  cat(sprintf("Error: input file: %s does not exist or cannot be read!\n\n", input_file))
  print_help(parser)
  quit()
}
if(is.null(opt$output_filename)){
  opt$output_filename = input_file
}

# read in alignments
alignments = read.table(input_file, stringsAsFactors = F, row.names=NULL, fill = T, header = F)[, c(1:12)]
# avoid inter overflow
alignments[, c(2:4, 7:12)] = apply(alignments[, c(2:4, 7:12)], 2, as.double)

# set column names
# PAF IS ZERO-BASED - CHECK HOW CODE WORKS
colnames(alignments)[1:12] = c("queryID","queryLen","queryStart","queryEnd","strand","refID","refLen","refStart","refEnd","numResidueMatches","lenAln","mapQ")

# caculate similarity
alignments$percentID = alignments$numResidueMatches / alignments$lenAln
if (opt$identity_lower) {
    alignments$percentID[which(alignments$percentID < opt$identity_lower)] <- opt$identity_lower
}

queryStartTemp = alignments$queryStart
# Flip starts, ends for negative strand alignments
alignments$queryStart[which(alignments$strand == "-")] = alignments$queryEnd[which(alignments$strand == "-")]
alignments$queryEnd[which(alignments$strand == "-")] = queryStartTemp[which(alignments$strand == "-")]
rm(queryStartTemp)

cat(paste0("\nNumber of alignments: ", nrow(alignments), "\n"))
cat(paste0("Number of query sequences: ", length(unique(alignments$queryID)), "\n"))

# sort by ref chromosome sizes, keep top X chromosomes OR keep specified IDs
if(is.null(opt$refIDs)){
  if (opt$sortbyID){
    refIDsToKeepOrdered = unique(sort(alignments$refID))
  }else{
    chromMax = tapply(alignments$refLen, alignments$refID, max)
    refIDsToKeepOrdered = names(sort(chromMax, decreasing = T))
  }
}else{
  refIDsToKeepOrdered = unlist(strsplit(opt$refIDs, ","))
  alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered),]
}

# filter queries by alignment length, for now include overlapping intervals
queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]),]
# filter alignment by length
alignments = alignments[which(alignments$lenAln > opt$min_align),]
# filter alignment by ref length
alignments = alignments[which(alignments$refLen > opt$min_ref_len),]
# re-filter queries by alignment length, for now include overlapping intervals
queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]),]

cat(paste0("\nAfter filtering... Number of alignments: ", nrow(alignments),"\n"))
cat(paste0("After filtering... Number of query sequences: ", length(unique(alignments$queryID)),"\n\n"))

if (nrow(alignments) == 0) {
  cat("Error: All alignments were discarded after filtering, no valid alignments remain.\n\n.")
  quit()
}

# filter refIDsToKeepOrdered if some refID are filtered
refIDsToKeepOrdered = refIDsToKeepOrdered[which(refIDsToKeepOrdered %in% alignments$refID)]

# sort df on ref
alignments$refID = factor(alignments$refID, levels = refIDsToKeepOrdered) # set order of refID
alignments = alignments[with(alignments,order(refID,refStart)),]
chromMax = tapply(alignments$refLen, alignments$refID, max)
# make new ref alignments for dot plot

alignments$refStart2 = alignments$refStart + sapply(as.character(alignments$refID), function(x) ifelse(x == names(chromMax)[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )
alignments$refEnd2 = alignments$refEnd + sapply(as.character(alignments$refID), function(x) ifelse(x == names(chromMax)[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )

## queryID sorting step 1/2
# sort levels of factor 'queryID' based on longest alignment
alignments$queryID = factor(alignments$queryID, levels=unique(as.character(alignments$queryID)))
queryMaxAlnIndex = tapply(alignments$lenAln, alignments$queryID, which.max, simplify = F)
alignments$queryID = factor(alignments$queryID, levels = unique(as.character(alignments$queryID))[order(mapply(
  function(x, i)
    alignments$refStart2[which(i == alignments$queryID)][x],
  queryMaxAlnIndex,
  names(queryMaxAlnIndex)
))])

## queryID sorting step 2/2
## sort levels of factor 'queryID' based on longest aggregrate alignmentst to refID's
# per query ID, get aggregrate alignment length to each refID 
queryLenAggPerRef = sapply((levels(alignments$queryID)), function(x) tapply(alignments$lenAln[which(alignments$queryID == x)], alignments$refID[which(alignments$queryID == x)], sum) )
if(length(levels(alignments$refID)) > 1){
  queryID_Ref = apply(queryLenAggPerRef, 2, function(x) rownames(queryLenAggPerRef)[which.max(x)])
} else {
  queryID_Ref = sapply(queryLenAggPerRef, function(x) names(queryLenAggPerRef)[which.max(x)])
}
# set order for queryID
alignments$queryID = factor(alignments$queryID, levels = (levels(alignments$queryID))[order(match(queryID_Ref, levels(alignments$refID)))])
queryMax = tapply(alignments$queryLen, alignments$queryID, max)

if(opt$flip){
  #  flip query starts stops to forward if most align are in reverse complement
  queryRevComp = tapply(alignments$queryEnd - alignments$queryStart, alignments$queryID, function(x) sum(x)) < 0
  queryRevComp = names(queryRevComp)[which(queryRevComp)]
  alignments$queryStart[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryStart[which(alignments$queryID %in% queryRevComp)] + 1
  alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] + 1
}
## make new query alignments for dot plot
alignments$queryStart2 = alignments$queryStart + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )
alignments$queryEnd2 = alignments$queryEnd + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )

# plot break points
if (opt$break_point) {
  alignments$break_col = rep(0, length(alignments$percentID));
}

options(warn = -1) # turn off warnings
gp = ggplot(alignments) + 
  theme_bw() + 
  theme(
    text = element_text(size = 12),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(angle = 15),
    axis.text.x = element_text(hjust = 1, angle = 45)
  ) +
  scale_color_distiller(palette = "Spectral") +
  labs(color = "Identity",
       title = paste0(paste0("Post-filtering number of alignments: ", nrow(alignments),"\t\t\t\t"),
                      paste0("minimum alignment length (-m): ", opt$min_align,"\n"),
                      paste0("Post-filtering number of queries: ", length(unique(alignments$queryID)),"\t\t\t\t\t\t\t\t"),
                      paste0("minimum query aggregate alignment length (-q): ", opt$min_query_aln)
       )
  )

# x
if (length(unique(alignments$refID)) == 1){
  reflen = unique(alignments$refLen)
  xbreaks = seq(0, reflen, reflen/10)
  if (reflen/10 > 1e9){
    xlables = paste(round(xbreaks/1e9), "GB")
  }else if (reflen/10 > 1e6) {
    xlables = paste(round(xbreaks/1e6), "MB")
  }else if(reflen/10 > 1e3){
    xlables = paste(round(xbreaks/1e3), "KB")
  }else{
    xlables = paste(round(xbreaks), "bp")
  }

  gp = gp + scale_x_continuous(expand = c(0, 0), limits = c(0, reflen + 0.1), breaks = xbreaks, labels = xlables) +
    xlab(unique(alignments$refID))
}else{
  gp = gp + 
    theme(panel.grid.major.x=element_blank()) +
    geom_vline(xintercept = cumsum(as.numeric(chromMax)), col="#ebebeb") + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, sum(as.numeric(chromMax)) + 0.1), 
      breaks = cumsum(as.numeric(chromMax)) - chromMax/2, labels = substr(levels(alignments$refID), start = 1, stop = 20)) + 
    xlab("Reference")
}

# y
if (length(unique(alignments$queryID)) == 1){
  queryLen = unique(alignments$queryLen)
  ybreaks = seq(0, queryLen, queryLen/10)
  if (queryLen/10 > 1e9){
    ylables = paste(round(ybreaks/1e9), "GB")
  }else if (queryLen/10 > 1e6) {
    ylables = paste(round(ybreaks/1e6), "MB")
  }else if(queryLen/10 > 1e3){
    ylables = paste(round(ybreaks/1e3), "KB")
  }else{
    ylables = paste(round(ybreaks), "bp")
  }

  gp = gp + scale_y_continuous(expand = c(0, 0), limits = c(0, queryLen + 0.1), breaks = ybreaks, labels = ylables) +
    ylab(unique(alignments$queryID))
}else{
    gp = gp + 
    theme(panel.grid.major.y=element_blank()) +
    geom_hline(yintercept = cumsum(as.numeric(queryMax)), col="#ebebeb") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, sum(as.numeric(queryMax)) + 0.1), 
      breaks = cumsum(as.numeric(queryMax)) - queryMax/2, labels = substr(levels(alignments$queryID), start = 1, stop = 20)) + 
    ylab("Query")
}

# plot co-line
gp = gp + geom_segment(aes(x = refStart2, xend = refEnd2, y = queryStart2, yend = queryEnd2, color = percentID))

# plot break points
if (opt$break_point) {
  gp = gp + geom_point(mapping = aes(x = refStart2, y = queryStart2, color = break_col), size = opt$plot_size/60, shape = 19) +
    geom_point(mapping = aes(x = refEnd2, y = queryEnd2, color = break_col), size = opt$plot_size/60, shape = 19)
}

# plot regions in ref bed
if (!is.null(opt$ref_bed_file)) {
  ref_bed = read.table(opt$ref_bed_file, header=F, stringsAsFactors=F)
  colnames(ref_bed)[1:3] = c("refID", "start", "end")
  ref_bed$color = generate_colors(nrow(ref_bed))

  ref_bed$xstart = ref_bed$start + sapply(ref_bed$refID, function(x) {
    ifelse(x == names(chromMax)[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax))-1])
  })
  ref_bed$xend = ref_bed$end + sapply(ref_bed$refID, function(x) {
    ifelse(x == names(chromMax)[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax))-1])
  })

  for (i in 1:nrow(ref_bed)) {
    gp = gp + geom_vline(xintercept = c(ref_bed$xstart[i], ref_bed$xend[i]), linetype="dashed", color=ref_bed$color[i])
  }
}

# plot regions in query bed
if (!is.null(opt$query_bed_file)) {
  query_bed = read.table(opt$query_bed_file, header=F, stringsAsFactors=F)
  colnames(query_bed)[1:3] = c("queryID", "start", "end")
  query_bed$color = generate_colors(nrow(query_bed))

  if (opt$flip) {
    rev_idx = which(query_bed$queryID %in% queryRevComp)
    if (length(rev_idx) > 0) {
      query_bed$start[rev_idx] = queryMax[query_bed$queryID[rev_idx]] - query_bed$start[rev_idx] + 1
      query_bed$end[rev_idx] = queryMax[query_bed$queryID[rev_idx]] - query_bed$end[rev_idx] + 1
    }
  }

  query_bed$ystart = query_bed$start + sapply(query_bed$queryID, function(x) {
    ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1])
  })
  query_bed$yend = query_bed$end + sapply(query_bed$queryID, function(x) {
    ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1])
  })

  for (i in 1:nrow(query_bed)) {
    gp = gp + geom_hline(yintercept = c(query_bed$ystart[i], query_bed$yend[i]),
                         linetype="twodash", color=query_bed$color[i])
  }
}

# save
ggsave(filename = paste0(opt$output_filename, ".svg"), width = opt$plot_size, height = opt$plot_size * 0.8, units = "in", dpi = 1200, limitsize = F)
options(warn=0) # turn on warnings
RS2
  chmod +x "$PAF2_R"

  if [[ -n "$REFERENCE" ]]; then
    PAF_OUT="${D5}/query_to_reference.paf"
    # log both commands to master log
    log "COMMAND: minimap2 -t \"$THREADS\" $MM2_PARAMS \"$REFERENCE\" \"$WORK_ASM\" > \"$PAF_OUT\""
    ( set -x; minimap2 -t "$THREADS" $MM2_PARAMS "$REFERENCE" "$WORK_ASM" > "$PAF_OUT" ) &>> "${PIPELINE_LOG}"
    log "COMMAND: \"$PAF2_R\" -f -b \"$PAF_OUT\""
    ( set -x; "$PAF2_R" -f -b "$PAF_OUT" ) &>> "${PIPELINE_LOG}"
  else
    log "Reference genome not provided (--reference). Skipping dot-plot generation."
  fi


  if [[ -s "${D5}/plotted_errors.png" && -s "${D5}/plotted_errors.svg" ]]; then
    write_marker "$D5" "step_5_done.txt" "ok"
  else
    write_marker "$D5" "step_5_done.txt" "failed"
  fi
fi

# --------------------------- finalization ---------------------------
log "Pipeline completed on $(date '+%Y-%m-%d %H:%M:%S')"
# total runtime
PIPELINE_END_TS="$(date +%s)"
PIPELINE_DUR_SEC=$((PIPELINE_END_TS - PIPELINE_START_TS))
printf -v __dur 'Total runtime: %02dh:%02dm:%02ds' $((PIPELINE_DUR_SEC/3600)) $(((PIPELINE_DUR_SEC%3600)/60)) $((PIPELINE_DUR_SEC%60))
log "$__dur"

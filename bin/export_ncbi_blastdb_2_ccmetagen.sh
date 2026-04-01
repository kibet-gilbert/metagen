#!/usr/bin/env bash
# Export sequences from a BLAST DB for CCMetagen:
#  - Export headers as >taxid|title (or >taxid|sci_name|title if you prefer)
#  - Title filtering (drop PREDICTED, vector/synthetic/patent/uncultured/environmental/metagenome)
#  - Length filter
#  - Deduplicate
#  - (Optional) low-complexity masking
#  - (Optional) compress final
# Logs include sequence counts at every step. Intermediates are removed by default.

set -Eeuo pipefail

# -------------------- User-tunable (env) --------------------
MIN_LEN="${MIN_LEN:-300}"                 # keep sequences >= this length
MASK_LOW_COMPLEX="${MASK_LOW_COMPLEX:-1}" # 1=run dustmasker; 0=skip
COMPRESS_FINAL="${COMPRESS_FINAL:-1}"     # 1=gzip final; 0=leave plain fasta
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}" # 1=keep tmp files; 0=delete on exit
THREADS="${THREADS:-4}"                   # threads for seqkit / dustmasker where applicable

# Header format for blastdbcmd:
# Default CCMetagen-friendly: >taxid|title\nsequence
# If you prefer to also include the scientific name: set HEADER_FMT=">%T|%S|%t\n%s"
HEADER_FMT="${HEADER_FMT:-">%T|%t\n%s"}"

# -------------------- Args & usage --------------------------
usage() {
  cat <<EOF
Usage: $0 <blast_db_prefix> [outdir]

Examples:
  $0 /export/data/bio/ncbi/blast/db/v5/nt                ./nt_out
  $0 /export/data/ilri/.../ncbi/core_nt/core_nt          ./ccm_out

Environment variables (override defaults):
  MIN_LEN=${MIN_LEN}
  MASK_LOW_COMPLEX=${MASK_LOW_COMPLEX}
  COMPRESS_FINAL=${COMPRESS_FINAL}
  KEEP_INTERMEDIATE=${KEEP_INTERMEDIATE}
  THREADS=${THREADS}
  HEADER_FMT='${HEADER_FMT}'
EOF
}

DB_PREFIX="${1:-}"
OUTDIR="${2:-ccm_out}"
[[ -z "${DB_PREFIX}" ]] && { usage; exit 1; }
mkdir -p "${OUTDIR}"

# -------------------- Logging helpers -----------------------
LOG="${OUTDIR}/run.log"
TMP="${OUTDIR}/.tmp"
mkdir -p "${TMP}"

ts()  { date +"%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*" | tee -a "${LOG}"; }
die() { echo "[$(ts)] ERROR: $*" | tee -a "${LOG}" >&2; exit 1; }

# -------------------- Tool checks ---------------------------
need() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }
need blastdbcmd
need seqkit
need awk; need sed; need grep
if [[ "${MASK_LOW_COMPLEX}" == "1" ]]; then need dustmasker; fi
if [[ "${COMPRESS_FINAL}" == "1" ]]; then need gzip; fi

# -------------------- File paths ----------------------------
RAW="${TMP}/01.raw.fasta"
CLEAN="${TMP}/02.clean_titles.fasta"
LEN="${TMP}/03.len${MIN_LEN}.fasta"
DEDUP="${TMP}/04.len${MIN_LEN}.dedup.fasta"
MASK="${TMP}/05.len${MIN_LEN}.dedup.masked.fasta"

FINAL="${OUTDIR}/ccmetagen.fasta"
[[ "${MASK_LOW_COMPLEX}" == "1" ]] && FINAL="${OUTDIR}/ccmetagen.masked.fasta"
FINAL_OUT="${FINAL}"

# -------------------- Utilities -----------------------------
count_fa() {    # count sequences in FASTA (plain or gz)
  local f="$1"
  if [[ -s "$f" ]]; then
    if [[ "$f" =~ \.gz$ ]]; then zgrep -c '^>' "$f" || echo 0
    else grep -c '^>' "$f" || echo 0
    fi
  else
    echo 0
  fi
}

cleanup() {
  if [[ "${KEEP_INTERMEDIATE}" != "1" ]]; then
    log "Cleaning intermediates in ${TMP}/"
    rm -rf "${TMP:?}/" || true
  else
    log "KEEP_INTERMEDIATE=1 — leaving intermediates in ${TMP}/"
  fi
}
trap cleanup EXIT

# Translate HEADER_FMT into a safe -outfmt (preserve \n)
# We’ll use ANSI-C quoting: $"...".
ansi_outfmt() {
  # Escape backslashes and quotes so the shell passes the string correctly
  printf "%s" "$HEADER_FMT" | sed -e 's/\\/\\\\/g' -e 's/"/\\"/g'
}

# -------------------- Start log -----------------------------
{
  echo "=============================================================="
  echo "Started        : $(ts)"
  echo "DB_PREFIX      : ${DB_PREFIX}"
  echo "OUTDIR         : ${OUTDIR}"
  echo "MIN_LEN        : ${MIN_LEN}"
  echo "MASK_LOW_COMPLEX: ${MASK_LOW_COMPLEX}"
  echo "COMPRESS_FINAL : ${COMPRESS_FINAL}"
  echo "KEEP_INTERMEDIATE: ${KEEP_INTERMEDIATE}"
  echo "THREADS        : ${THREADS}"
  echo "HEADER_FMT     : ${HEADER_FMT}"
  echo "=============================================================="
} | tee -a "${LOG}"

# -------------------- STEP 1: Export ------------------------
log "Step 1/5: Exporting from BLAST DB with headers ${HEADER_FMT}"
OFMT=$(ansi_outfmt)
t0=$(date +%s)

blastdbcmd -db "${DB_PREFIX}" \
  -entry all \
  -outfmt $"$OFMT" \
  > "${RAW}"

t1=$(date +%s)
log "  Export time: $((t1 - t0))s"
log "  Exported sequences: $(count_fa "${RAW}")"
log "  RAW: ${RAW}"

# -------------------- STEP 2: Title filtering ---------------
log "Step 2/5: Filtering titles (drop PREDICTED / vector / synthetic construct / patent / uncultured / environmental sample / metagenome)"
export LC_ALL=C

# Keep only records whose header DOES NOT match the excluded patterns:
#   - '^>[^|]+|PREDICTED:'  -> >taxid|PREDICTED:...
#   - titles containing: vector|cloning|synthetic construct|patent|uncultured|environmental sample|metagenome
grep -Eiv '^[^>]*>[^|]+\|PREDICTED:' "${RAW}" \
  | grep -Eiv 'vector|cloning|synthetic construct|patent|uncultured|environmental sample|metagenome' \
  | awk 'BEGIN{RS=">"; ORS=""} NR>1 {print ">"$0}' \
  > "${CLEAN}"

log "  After title filter: $(count_fa "${CLEAN}")"
log "  CLEAN: ${CLEAN}"

# -------------------- STEP 3: Length filter -----------------
log "Step 3/5: Length filter (>= ${MIN_LEN} nt)"
seqkit seq -j "${THREADS}" -m "${MIN_LEN}" "${CLEAN}" > "${LEN}"
log "  After length filter: $(count_fa "${LEN}")"
log "  LEN: ${LEN}"

# -------------------- STEP 4: Deduplicate -------------------
log "Step 4/5: Deduplicate identical sequences (seqkit rmdup -s -i)"
seqkit rmdup -j "${THREADS}" -s -i "${LEN}" > "${DEDUP}"
log "  After dedup: $(count_fa "${DEDUP}")"
log "  DEDUP: ${DEDUP}"

# -------------------- STEP 5: Masking (optional) ------------
if [[ "${MASK_LOW_COMPLEX}" == "1" ]]; then
  log "Step 5/5: Low-complexity masking (dustmasker -> hard-mask to N)"
  dustmasker -in "${DEDUP}" -outfmt fasta \
    | sed 's/[a-z]/N/g' \
    > "${MASK}"
  log "  After masking: $(count_fa "${MASK}")"
  cp -f "${MASK}" "${FINAL}"
else
  log "Step 5/5: Skipping masking (MASK_LOW_COMPLEX=0)"
  cp -f "${DEDUP}" "${FINAL}"
fi

# -------------------- Finalize ------------------------------
if [[ "${COMPRESS_FINAL}" == "1" ]]; then
  log "Compressing final FASTA (gzip)"
  gzip -f "${FINAL}"
  FINAL_OUT="${FINAL}.gz"
fi

log "FINAL: ${FINAL_OUT}"
log "Sequences in final: $(count_fa "${FINAL_OUT}")"
log "Done: $(ts)"

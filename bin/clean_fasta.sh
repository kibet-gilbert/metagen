#!/usr/bin/env bash
# =============================================================================
# clean_fasta.sh   —   MODULE 2 of 2
#
# Clean a raw FASTA file (from export_blastdb.sh OR download_all_refseq.sh)
# for use as a CCMetagen / KMA database:
#
#   Step 1 — Title filtering  (drop PREDICTED, vectors, synthetic, uncultured...)
#   Step 2 — Length filter    (keep sequences >= MIN_LEN)
#   Step 3 — Deduplication    (seqkit rmdup by sequence)
#   Step 4 — Low-complexity masking (dustmasker → hard-mask to N)  [optional]
#   Step 5 — Compress         (pigz/gzip)                          [optional]
#
# The script accepts plain or gzip-compressed input.
# Intermediate files are removed unless KEEP_INTERMEDIATE=1.
#
# Usage:
#   clean_fasta.sh <input.fasta[.gz]> [outdir]
#
# Examples:
#   # After export_blastdb.sh:
#   clean_fasta.sh ./nt_export/nt.raw.fasta.gz  ./nt_clean
#
#   # After download_all_refseq.sh (merged RefSeq):
#   clean_fasta.sh ./refseq/refseq_ccmetagen.fna.gz  ./refseq_clean
#
# Environment overrides:
#   MIN_LEN=300              minimum sequence length (default: 300)
#   MASK_LOW_COMPLEX=1       1=run dustmasker (default); 0=skip
#   DEDUP_SEQUENCES=0        1=Filter Duplicate Sequences (default); 0=skip
#   COMPRESS_FINAL=1         1=gzip output (default); 0=plain
#   KEEP_INTERMEDIATE=0      1=keep tmp files; 0=delete (default)
#   THREADS=4                threads for seqkit / pigz
# =============================================================================
set -Eeuo pipefail

on_error() {
  local exit_code=$?
  local line_no=$1
  echo "Error: exit code ${exit_code} at line ${line_no}" >&2
}

trap 'on_error $LINENO' ERR

# ── Defaults ──────────────────────────────────────────────────────────────────
MIN_LEN="${MIN_LEN:-300}"
MASK_LOW_COMPLEX="${MASK_LOW_COMPLEX:-0}"
DEDUP_SEQUENCES="${DEDUP_SEQUENCES:-0}"
COMPRESS_FINAL="${COMPRESS_FINAL:-1}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"
THREADS="${THREADS:-4}"

usage() {
  cat <<EOF
Usage: $0 <input.fasta[.gz]> [outdir]

  input   : raw FASTA (plain or .gz) — from export_blastdb.sh or download_all_refseq.sh
  outdir  : output directory (default: ./fasta_clean)

Environment:
  MIN_LEN=${MIN_LEN}
  MASK_LOW_COMPLEX=${MASK_LOW_COMPLEX}   (1=dustmasker, 0=skip)
  DEDUP_SEQUENCES=${DEDUP_SEQUENCES}     (1=sedkit rmdup, 0=skip)
  COMPRESS_FINAL=${COMPRESS_FINAL}       (1=gzip, 0=plain)
  KEEP_INTERMEDIATE=${KEEP_INTERMEDIATE} (1=keep tmp, 0=delete)
  THREADS=${THREADS}
EOF
}

INPUT="${1:-}"
OUTDIR="${2:-fasta_clean}"
[[ -z "$INPUT" ]] && { usage; exit 1; }
[[ ! -f "$INPUT" ]] && { echo "ERROR: Input not found: $INPUT"; exit 1; }
mkdir -p "$OUTDIR"

# ── Logging ───────────────────────────────────────────────────────────────────
LOG="${OUTDIR}/clean_fasta.log"
TMP="${OUTDIR}/.tmp"
mkdir -p "$TMP"

ts()  { date +"%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*" | tee -a "$LOG"; }
die() { echo "[$(ts)] ERROR: $*" | tee -a "$LOG" >&2; exit 1; }

# ── Tool checks ───────────────────────────────────────────────────────────────
for t in seqkit awk sed grep; do
  command -v "$t" >/dev/null 2>&1 || die "Missing: $t"
done
[[ "$MASK_LOW_COMPLEX" == "1" ]] && \
  command -v dustmasker >/dev/null 2>&1 || { log "[WARN] dustmasker not found — disabling masking"; MASK_LOW_COMPLEX=0; }
[[ "$COMPRESS_FINAL" == "1" ]] && \
  { PZIP=$(command -v pigz 2>/dev/null || command -v gzip 2>/dev/null || die "pigz/gzip not found"); }

# ── Helper: count sequences ───────────────────────────────────────────────────
count_fa() {
  local f="$1"
  if [[ "$f" =~ \.gz$ ]]; then zgrep -c '^>' "$f" 2>/dev/null || echo 0
  else grep -c '^>' "$f" 2>/dev/null || echo 0
  fi
}

# ── Helper: cat plain or gz ───────────────────────────────────────────────────
cat_fa() {
  if [[ "$1" =~ \.gz$ ]]; then zcat "$1"; else cat "$1"; fi
}

# ── Intermediates ─────────────────────────────────────────────────────────────
STEP1="${TMP}/01.title_filtered.fasta"
STEP2="${TMP}/02.len${MIN_LEN}.fasta"
STEP3="${TMP}/03.dedup.fasta"
STEP4="${TMP}/04.masked.fasta"

BASENAME=$(basename "${INPUT%.gz}")
BASENAME="${BASENAME%.fasta}"
BASENAME="${BASENAME%.fna}"
FINAL="${OUTDIR}/${BASENAME}.ccmetagen.fasta"
[[ "$MASK_LOW_COMPLEX" == "1" ]] && FINAL="${OUTDIR}/${BASENAME}.ccmetagen.masked.fasta"
[[ "$COMPRESS_FINAL"   == "1" ]] && FINAL="${FINAL}.gz"

# ── Cleanup on exit ───────────────────────────────────────────────────────────
cleanup() {
  if [[ "$KEEP_INTERMEDIATE" != "1" ]]; then
    log "Cleaning intermediates in $TMP/"
    rm -rf "$TMP" || true
  else
    log "KEEP_INTERMEDIATE=1 — leaving $TMP/ intact"
  fi
}
trap cleanup EXIT

{
  echo "======================================================================"
  echo "clean_fasta.sh — MODULE 2"
  echo "Started    : $(ts)"
  echo "INPUT      : $INPUT"
  echo "OUTDIR     : $OUTDIR"
  echo "MIN_LEN    : $MIN_LEN"
  echo "MASK       : $MASK_LOW_COMPLEX"
  echo "DEDUPLICATE: $DEDUP_SEQUENCES"
  echo "COMPRESS   : $COMPRESS_FINAL"
  echo "THREADS    : $THREADS"
  echo "======================================================================"
} | tee -a "$LOG"

n_INPUT=$(count_fa "$INPUT")
log "Input sequences: ${n_INPUT}"

# ── STEP 1: Title filtering ───────────────────────────────────────────────────
log "Step 1/4: Title filtering"
log "  Dropping: vector | synthetic construct | cloning | patent | metagenome"

export LC_ALL=C

# Strategy: convert FASTA to per-record awk blocks, filter on header, then
# reconstruct. Works on both plain and gzipped input via cat_fa.
cat_fa "$INPUT" \
  | awk '
      /^>/ { header = $0; seq = ""; next }
      { seq = seq $0 }
      END { }          # handled below
    ' 2>/dev/null || true   # awk streaming alternative below

# Simpler streaming approach that handles multi-line FASTA:
cat_fa "$INPUT" \
  | awk 'BEGIN{RS=">"; ORS=""} NR>1 {print ">"$0}' \
  | grep -Eiv 'vector|synthetic construct|cloning|patent|metagenome' \
  | awk 'BEGIN{RS=">"; ORS=""} NR>1 {print ">"$0}' \
  > "$STEP1"

n_STEP1=$(count_fa "$STEP1")
log "  After title filter: ${n_STEP1}"

# ── STEP 2: Length filter ─────────────────────────────────────────────────────
log "Step 2/4: Length filter (>= ${MIN_LEN} nt)"
#[WARN] you may switch on flag -g/--remove-gaps to remove spaces
seqkit seq -j "$THREADS" -m "$MIN_LEN" "$STEP1" > "$STEP2"
n_STEP2=$(count_fa "$STEP2")
log "  After length filter: ${n_STEP2}"

# ── STEP 3: Deduplicate ───────────────────────────────────────────────────────
log "Step 3/4: Identify identical sequences (seqkit rmdup -s -i)"
# Run once
seqkit rmdup -j "$THREADS" -s -i "$STEP2" \
  --id-regexp '^([^|]+\|[A-Za-z.]+(?:\s+[A-Za-z]+)?)' \
  -d "${OUTDIR}/duplicates.txt" \
  -D "${OUTDIR}/duplicate_seqids.txt" \
  > "$STEP3"

# Precompute counts (avoid repeated expensive calls)
n_STEP3=$(count_fa "$STEP3")
n_duplicates=$(( n_STEP2 - n_STEP3 ))

if [[ "$DEDUP_SEQUENCES" == "1" ]]; then
  log "  Deduplicated identical sequences"
  log "  After dedup: ${n_STEP3} (removed ${n_duplicates} duplicates)"
  STEP4_INPUT=$STEP3
  n_STEP4_INPUT=${n_STEP3}
else
  log "  Deduplication skipped (reporting only)"
  log "  Found ${n_duplicates} duplicates)"
  STEP4_INPUT=$STEP2
  n_STEP4_INPUT=${n_STEP2}
fi
log "  Duplicate details: ${OUTDIR}/duplicates.txt"

# ── STEP 4: Low-complexity masking ───────────────────────────────────────────
if [[ "$MASK_LOW_COMPLEX" == "1" ]]; then
  log "Step 4/4: Low-complexity masking (dustmasker → hard-mask lowercase → N)"
  dustmasker -in "$STEP4_INPUT" -outfmt fasta \
    | awk '
        /^>/ {print; next}
        {gsub(/[a-z]/,"N"); print}
      ' \
    > "$STEP4"
  n_STEP4=$(count_fa "$STEP4")
  log "  Masking complete: ${n_STEP4}"
  PRE_COMPRESS="$STEP4"
else
  log "Step 4/4: Masking skipped (MASK_LOW_COMPLEX=0)"
  PRE_COMPRESS="$STEP4_INPUT"
fi

# ── Compress and write final ──────────────────────────────────────────────────
log "Writing final output: $FINAL"
if [[ "$COMPRESS_FINAL" == "1" ]]; then
  $PZIP -p "$THREADS" -c "$PRE_COMPRESS" > "$FINAL"
else
  cp -f "$PRE_COMPRESS" "$FINAL"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
n_final=${n_STEP4}
sz_final=$(du -sh "$FINAL" | cut -f1)

{
  echo "======================================================================"
  echo "clean_fasta.sh SUMMARY"
  echo "Input sequences     : ${n_INPUT}"
  echo "After title filter  : ${n_STEP1}"
  echo "After length (${MIN_LEN}): ${n_STEP2}"
  echo "After dedup         : ${n_STEP4_INPUT}"
  echo "Duplicate sequences : ${n_duplicates}"
  echo "Final output        : $n_final sequences  ($sz_final)"
  echo "Output file         : $FINAL"
  echo "Finished            : $(ts)"
  echo "======================================================================"
} | tee -a "$LOG"

log ""
log "Next step — build KMA index:"
log "  kma index -i $FINAL -o <db_prefix>"
log ""
log "Done: $(ts)"

#!/usr/bin/env bash
# =============================================================================
# export_blastdb.sh   —   MODULE 1 of 2
#
# Export raw sequences from an NCBI BLAST database (nt, core_nt, etc.)
# with CCMetagen-compatible headers: >taxid|title
#
# Output: a single gzip-compressed FASTA ready to pipe into clean_fasta.sh
#
# Usage:
#   export_blastdb.sh <blast_db_prefix> [outdir]
#
# Examples:
#   export_blastdb.sh /export/data/ncbi/blast/db/v5/nt       ./nt_export
#   export_blastdb.sh /export/data/ncbi/blast/db/v5/core_nt  ./core_nt_export
#
# Environment overrides:
#   HEADER_FMT   blastdbcmd -outfmt format string
#                Default: >%T|%t\n%s   (taxid|title)
#   THREADS      threads for pigz compression  (default: 4)
#   COMPRESS     1 = gzip output (default); 0 = plain FASTA
# =============================================================================
set -Eeuo pipefail

on_error() {
  local exit_code=$?
  local line_no=$1
  echo "Error: exit code ${exit_code} at line ${line_no}" >&2
}

trap 'on_error $LINENO' ERR
# =============================================================================

HEADER_FMT="${HEADER_FMT:->%T|%t\n%s}"
THREADS="${THREADS:-4}"
COMPRESS="${COMPRESS:-1}"

usage() {
  cat <<EOF
Usage: $0 <blast_db_prefix> [outdir]

  blast_db_prefix : path to BLAST DB without extension (e.g. /data/ncbi/nt)
  outdir          : output directory (default: ./blast_export)

Environment:
  HEADER_FMT='${HEADER_FMT}'   (blastdbcmd -outfmt string)
  THREADS=${THREADS}
  COMPRESS=${COMPRESS}           (1=gzip, 0=plain)
EOF
}

DB_PREFIX="${1:-}"
OUTDIR="${2:-blast_export}"
[[ -z "$DB_PREFIX" ]] && { usage; exit 1; }
mkdir -p "$OUTDIR"

# ── Logging ───────────────────────────────────────────────────────────────────
LOG="${OUTDIR}/export_blastdb.log"
ts()  { date +"%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*" | tee -a "$LOG"; }
die() { echo "[$(ts)] ERROR: $*" | tee -a "$LOG" >&2; exit 1; }

# ── Tool checks ───────────────────────────────────────────────────────────────
command -v blastdbcmd >/dev/null 2>&1 || die "blastdbcmd not found (load BLAST module)"
[[ "$COMPRESS" == "1" ]] && \
  { command -v pigz >/dev/null 2>&1 || command -v gzip >/dev/null 2>&1 || die "pigz/gzip not found"; }

# ── Output paths ──────────────────────────────────────────────────────────────
DB_NAME=$(basename "$DB_PREFIX")
RAW_FASTA="${OUTDIR}/${DB_NAME}.raw.fasta"
[[ "$COMPRESS" == "1" ]] && RAW_FASTA="${RAW_FASTA}.gz"

{
  echo "======================================================================"
  echo "export_blastdb.sh — MODULE 1"
  echo "Started    : $(ts)"
  echo "DB_PREFIX  : $DB_PREFIX"
  echo "OUTDIR     : $OUTDIR"
  echo "HEADER_FMT : $HEADER_FMT"
  echo "COMPRESS   : $COMPRESS"
  echo "THREADS    : $THREADS"
  echo "======================================================================"
} | tee -a "$LOG"

# ── Verify DB exists ──────────────────────────────────────────────────────────
blastdbcmd -db "$DB_PREFIX" -info >/dev/null 2>&1 \
  || die "Cannot open BLAST database: $DB_PREFIX"

log "DB info:"
blastdbcmd -db "$DB_PREFIX" -info 2>&1 | tee -a "$LOG"

# ── Export ────────────────────────────────────────────────────────────────────
log "Exporting sequences with header format: $HEADER_FMT"
t0=$(date +%s)

if [[ "$COMPRESS" == "1" ]]; then
  PZIP=$(command -v pigz 2>/dev/null || echo "gzip")
  blastdbcmd -db "$DB_PREFIX" \
    -entry all \
    -outfmt $"$HEADER_FMT" \
    | $PZIP -p "$THREADS" \
    > "$RAW_FASTA"
else
  blastdbcmd -db "$DB_PREFIX" \
    -entry all \
    -outfmt $"$HEADER_FMT" \
    > "$RAW_FASTA"
fi

t1=$(date +%s)
log "Export complete: $((t1 - t0))s"

# ── Count sequences ───────────────────────────────────────────────────────────
if [[ "$COMPRESS" == "1" ]]; then
  n_seq=$(zgrep -c '^>' "$RAW_FASTA" || echo 0)
else
  n_seq=$(grep  -c '^>' "$RAW_FASTA" || echo 0)
fi
n_size=$(du -sh "$RAW_FASTA" | cut -f1)

log "Sequences exported : $n_seq"
log "File size          : $n_size"
log "Output             : $RAW_FASTA"
log ""
log "Next step:"
log "  Pass this file to clean_fasta.sh (Module 2):"
log "  clean_fasta.sh $RAW_FASTA $OUTDIR/cleaned"
log ""
log "Done: $(ts)"

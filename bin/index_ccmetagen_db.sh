#!/bin/bash
# =============================================================================
# index_ccmetagen_db.sh
#
# Build a KMA index from a pre-cleaned, taxonomy-annotated CCMetagen FASTA
# and create a compressed archive ready for local reuse or Zenodo upload.
#
# Usage:
# sbatch --job-name=kma_index_core_nt \
#        --partition=batch --nodelist 'compute05' \
#        --ntasks=1 --cpus-per-task=4 --mem=200G --time=12:00:00 \
#        --output=logs/kma_index_core_nt_%j.out \
#        --error=logs/kma_index_core_nt_%j.err \
#        --wrap="srun bash ~/bioinformatics/github/metagenomics/scripts/metagen/bin/index_ccmetagen_db.sh \
#        --fasta  core_nt/core_nt.taxid_desc.ccmetagen.fasta.gz \
#        --db-id  corent_kma \
#        --outdir core_nt/ --kma-args '-NI -Sparse TG -ME -CS 100000000'"
#
# The script will create:
#   kma_indexed_dbs/
#   └── core_nt/corent_kma
#       ├── core_nt.comp.b
#       ├── core_nt.length.b
#       ├── core_nt.name
#       ├── core_nt.seq.b
#       └── ...
#   core_nt/corent_kma.tar.gz   ← upload this to Zenodo
#
# Requirements: kma >= 1.4.x  (bioconda::kma)
# =============================================================================

set -euo pipefail
# IFS=$'\n\t'

# ── Defaults ──────────────────────────────────────────────────────────────────
FASTA=""
DB_ID=""
OUTDIR="kma_indexed_dbs"
KMA_ARGS_STR=""          # raw string from --kma-args
declare -a EXTRA_KMA_ARGS=()
# EXTRA_KMA_ARGS="-NI -Sparse TG -ME -CS 100000000"
LOGFILE=""

# ── Argument parsing ──────────────────────────────────────────────────────────
print_usage() {
    cat <<EOF

Usage:
  $(basename "$0") --fasta <file.fasta.gz> --db-id <name> [--outdir <dir>] [-- <kma index args>]

Options:
  --fasta   <file>   Input FASTA file (gzip accepted). Required.
  --db-id   <name>   Short identifier used as index prefix and archive name.
                     E.g. 'core_nt' or 'nt'. Required.
  --outdir  <dir>    Parent output directory [default: kma_indexed_dbs]
  --kma-args "<flags>"  Extra flags passed verbatim to kma_index as a quoted
                     string. Do NOT use -- for this (srun intercepts it).
                     E.g. --kma-args "-ME -CS 100000000 -Sparse TG"

Examples:
  # Index the core_nt CCMetagen FASTA
  ./index_ccmetagen_db.sh \\
      --fasta ccmetagen_db/core_nt/core_nt.taxid_desc.ccmetagen.fasta.gz \\
      --db-id core_nt

  # Index the full nt FASTA with a custom batch size
  ./index_ccmetagen_db.sh \\
      --fasta ccmetagen_db/nt/blast_nt.taxid_desc.ccmetagen.fasta.gz \\
      --db-id nt \\
      --outdir /var/scratch/gkibet/kma_dbs

  # Pass extra kma index flags after --
  ./index_ccmetagen_db.sh \\
      --fasta cleaned.fasta.gz --db-id my_db --kma-args "-k 16 -NI -Sparse TG -ME -CS 100000000"

EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta)   FASTA="$2";  shift 2 ;;
        --db-id)   DB_ID="$2";  shift 2 ;;
        --outdir)  OUTDIR="$2"; shift 2 ;;
        -h|--help) print_usage; exit 0  ;;
        --kma-args) KMA_ARGS_STR="$2"; shift 2 ;;
        *)         echo "[ERROR] Unknown option: $1"; print_usage; exit 1 ;;
    esac
done

# ── Load module ──────────────────────────────────────────────────────────────
module load kma/1.4.21

# ── Validation ────────────────────────────────────────────────────────────────
# Split --kma-args string into array elements (safe: no srun/-- mangling)
if [[ -n "$KMA_ARGS_STR" ]]; then
    read -ra EXTRA_KMA_ARGS <<< "$KMA_ARGS_STR"
fi

[[ -z "$FASTA" ]]  && { echo "[ERROR] --fasta is required.";  print_usage; exit 1; }
[[ -z "$DB_ID" ]]  && { echo "[ERROR] --db-id is required.";  print_usage; exit 1; }
[[ ! -f "$FASTA" ]] && { echo "[ERROR] FASTA file not found: $FASTA"; exit 1; }
command -v kma &>/dev/null || { echo "[ERROR] kma not found in PATH. Install with: conda install -c bioconda kma"; exit 1; }

# ── Setup ─────────────────────────────────────────────────────────────────────
DB_DIR="${OUTDIR}/${DB_ID}"
DB_PREFIX="${DB_DIR}/${DB_ID}"
ARCHIVE="${OUTDIR}/${DB_ID}.tar.gz"
LOGFILE="${OUTDIR}/${DB_ID}_index.log"
START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
FASTA_ABS=$(realpath "$FASTA")

mkdir -p "$DB_DIR"

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOGFILE"; }

log "======================================================================="
log " CCMetagen KMA Index Builder"
log "======================================================================="
log " FASTA input   : ${FASTA_ABS}"
log " DB identifier : ${DB_ID}"
log " Output prefix : ${DB_PREFIX}"
log " Archive output: ${ARCHIVE}"
log " KMA version   : $(kma -v 2>&1 | head -1)"
log " Start time    : ${START_TIME}"
log " Extra kma args: ${EXTRA_KMA_ARGS[*]:-(none)}"
log "======================================================================="

# ── FASTA stats (quick sanity check before committing hours of indexing) ──────
# BUG FIX: sequence counting via zcat | grep -c on a 305G gzip took 2.5 hours
# before indexing even started. Skip the count for files > SIZE_LIMIT_GB GB
# and report the file size only — still a useful sanity check, costs seconds.
SIZE_LIMIT_GB=10
FASTA_SIZE=$(du -sh "$FASTA" | cut -f1)
FASTA_BYTES=$(du -b  "$FASTA" | cut -f1)
SIZE_LIMIT_BYTES=$(( SIZE_LIMIT_GB * 1024 * 1024 * 1024 ))
 
if (( FASTA_BYTES > SIZE_LIMIT_BYTES )); then
    log "FASTA size     : ${FASTA_SIZE}  (> ${SIZE_LIMIT_GB}G — skipping sequence count to save time)"
    N_SEQS="(skipped for large file)"
else
    log "Counting input sequences (file < ${SIZE_LIMIT_GB}G, this is fast)..."
    if [[ "$FASTA" == *.gz ]]; then
        N_SEQS=$(zcat "$FASTA" | grep -c '^>' || true)
    else
        N_SEQS=$(grep -c '^>' "$FASTA" || true)
    fi
    log "  Sequences    : ${N_SEQS}"
    log "  File size    : ${FASTA_SIZE}"
fi

# ── KMA index ─────────────────────────────────────────────────────────────────
# Notes on kma index:
#   • Single-threaded — no -t flag (unlike kma align).
#   • Reads gzip input natively — no need to decompress first.
#   • Output: <prefix>.comp.b  .length.b  .name  .seq.b  (and .index.b for some versions)
#   • For CCMetagen, sequence headers MUST follow the format:
#       >taxid|<ncbi_taxid>|<description>
#     The ccmetagen.fasta files in your ccmetagen_db/ already have this format.
# command:
KMA_CMD=(
    kma_index
    -i "${FASTA_ABS}"
    -o "${DB_PREFIX}"
)
(( ${#EXTRA_KMA_ARGS[@]} > 0 )) && KMA_CMD+=("${EXTRA_KMA_ARGS[@]}")

log "Running kma index..."
log "  Command: ${KMA_CMD[*]}"
log "  Args breakdown:"
for arg in "${KMA_CMD[@]}"; do
    log "    [${arg}]"
done
 
"${KMA_CMD[@]}" 2>&1 | tee -a "$LOGFILE"
 
log "KMA index complete."
log "Index files created:"
ls -lh "${DB_DIR}/"* | tee -a "$LOGFILE"

# ── Checksum ──────────────────────────────────────────────────────────────────
log "Generating MD5 checksums..."
cd "$DB_DIR"
md5sum ./* > "${DB_ID}.md5"
cd - > /dev/null
log "Checksums written to ${DB_DIR}/${DB_ID}.md5"

# ── Metadata file (important for Zenodo) ─────────────────────────────────────
cat > "${DB_DIR}/README.txt" <<EOF
CCMetagen KMA Database — ${DB_ID}
==================================

Built   : ${START_TIME}
FASTA   : $(basename "${FASTA_ABS}")
Seqs    : ${N_SEQS}
KMA     : $(kma -v 2>&1 | head -1)
Builder : $(hostname)

Index prefix
  ${DB_ID}.*

Required files
  ${DB_ID}.comp.b     Compressed k-mer composition
  ${DB_ID}.length.b   Sequence lengths
  ${DB_ID}.name       Sequence names (taxid|desc format)
  ${DB_ID}.seq.b      Encoded sequences

Usage in CCMetagen / Nextflow
  ccmetagen -i sample.fq -o results -db /path/to/${DB_ID}/${DB_ID}

  # or via kma directly:
  kma -i sample_R1.fq sample_R2.fq \\
      -o results/${DB_ID} \\
      -t_db /path/to/${DB_ID}/${DB_ID} \\
      -1t1 -mem_mode -and -apm f

Checksums
  See ${DB_ID}.md5

Reference
  Dam P, Kania SA, Sørensen SJ et al. (2021). CCMetagen: a metagenomics
  tool based on KMA for fast and accurate identification of species in
  metagenomic datasets. Genome Biology 22, 122.
  https://doi.org/10.1186/s13059-021-02312-5

EOF

# ── Archive ───────────────────────────────────────────────────────────────────
log "Creating archive: ${ARCHIVE}"
tar -czf "$ARCHIVE" \
    -C "$OUTDIR" \
    "$(basename "$DB_DIR")"

ARCHIVE_SIZE=$(du -sh "$ARCHIVE" | cut -f1)
log "Archive created: ${ARCHIVE}  (${ARCHIVE_SIZE})"

# Final MD5 of the archive (needed for Zenodo upload verification)
md5sum "$ARCHIVE" > "${ARCHIVE}.md5"
log "Archive MD5: $(cat "${ARCHIVE}.md5")"

END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
log "======================================================================="
log " Done. End time: ${END_TIME}"
log " Archive to upload: ${ARCHIVE}"
log " Archive MD5      : ${ARCHIVE}.md5"
log "======================================================================="

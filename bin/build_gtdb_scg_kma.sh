#!/usr/bin/env bash
# =============================================================================
# build_gtdb_scg_kma.sh
#
# Extract, clean, and prepare GTDB bac120 single-copy marker gene (SCG)
# nucleotide sequences for KMA indexing and Genome Equivalent (GEQ)
# normalisation of ARG abundance.
#
# Usage:
#   build_gtdb_scg_kma.sh [gtdb_dir] [outdir] [reps|all] [threads]
#
# Environment overrides (all optional):
#   GTDB_DIR, OUTDIR, USE_REPS, THREADS
#   KEEP_EXTRACTED=1   keep per-gene extracted files
#   SKIP_CLEAN=1       skip seqkit length filter + dedup
#   MIN_LEN=100        minimum nucleotide length (default 100)
#   DEDUPLICATE_SEQUENCES=1	Deduplicate sequences(default: 0=skip)
#   SCG_FASTA_BUILDER=/path/to/build_gtdb_scg_fasta.py
#
# References:
#   Parks et al. (2022) GTDB: ongoing census. NAR 50:D785-D794.
#   Clausen et al. (2018) KMA. BMC Bioinformatics 19:307.
#   Nayfach & Bhatt (2016) Microbial GEQ. Nat Methods 13:290-292.
# =============================================================================
set -Eeuo pipefail

on_error() {
  local exit_code=$?
  local line_no=$1
  echo "Error: exit code ${exit_code} at line ${line_no}" >&2
}

trap 'on_error $LINENO' ERR
# =============================================================================

GTDB_DIR="${1:-${GTDB_DIR:-./gtdb_input}}"
OUTDIR="${2:-${OUTDIR:-./gtdb_scg_kma}}"
USE_REPS="${3:-${USE_REPS:-reps}}"
THREADS="${4:-${THREADS:-8}}"
KEEP_EXTRACTED="${KEEP_EXTRACTED:-0}"
SKIP_CLEAN="${SKIP_CLEAN:-0}"
MIN_LEN="${MIN_LEN:-100}"
DEDUPLICATE_SEQUENCES="${DEDUPLICATE_SEQUENCES:-0}"
SCG_FASTA_BUILDER="${SCG_FASTA_BUILDER:-./build_gtdb_scg_fasta.py}"

mkdir -p "$OUTDIR"

LOG="${OUTDIR}/build_gtdb_scg_kma.log"
ts()  { date +"%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*" | tee -a "$LOG"; }
die() { echo "[$(ts)] ERROR: $*" | tee -a "$LOG" >&2; exit 1; }

{
  echo "======================================================================"
  echo "build_gtdb_scg_kma.sh  started $(ts)"
  echo "GTDB_DIR  : $GTDB_DIR"
  echo "OUTDIR    : $OUTDIR"
  echo "USE_REPS  : $USE_REPS  (reps=species-reps only | all=all GTDB genomes)"
  echo "THREADS   : $THREADS"
  echo "MIN_LEN   : $MIN_LEN nt"
  echo "======================================================================"
} | tee -a "$LOG"

for t in tar awk sed python3 pigz; do
  command -v "$t" >/dev/null 2>&1 || die "Missing: $t"
done
HAS_SEQKIT=1; command -v seqkit >/dev/null 2>&1 || { log "[WARN] seqkit absent — skipping clean"; HAS_SEQKIT=0; SKIP_CLEAN=1; }
HAS_KMA=1;    command -v kma    >/dev/null 2>&1 || { log "[INFO] kma absent — skipping index"; HAS_KMA=0; }
[[ -x "$SCG_FASTA_BUILDER" ]] || die "Missing or non-executable: $SCG_FASTA_BUILDER"

# ── Identify archive ─────────────────────────────────────────────────────────
if [[ "$USE_REPS" == "reps" ]]; then
  ARCHIVE=$(ls "${GTDB_DIR}"/bac120_marker_genes_reps_r*.tar.gz 2>/dev/null | head -1 || true)
  LABEL="reps"
else
  ARCHIVE=$(ls "${GTDB_DIR}"/bac120_marker_genes_all_r*.tar.gz  2>/dev/null | head -1 || true)
  LABEL="all"
fi
[[ -z "$ARCHIVE" || ! -f "$ARCHIVE" ]] && \
  die "Cannot find bac120_marker_genes_${USE_REPS}_r*.tar.gz in $GTDB_DIR"
log "Archive : $ARCHIVE ($(du -sh "$ARCHIVE" | cut -f1))"

MARKER_INFO=$(ls  "${GTDB_DIR}"/bac120_msa_marker_info_r*.tsv 2>/dev/null | head -1 || true)
TAX_FILE=$(   ls  "${GTDB_DIR}"/bac120_taxonomy_r*.tsv*       2>/dev/null | head -1 || true)
[[ -f "$MARKER_INFO" ]] && { cp -f "$MARKER_INFO" "${OUTDIR}/marker_info.tsv"; log "Marker info: $MARKER_INFO"; }
if [[ -f "$TAX_FILE" ]]; then
  [[ "$TAX_FILE" =~ \.gz$ ]] && zcat "$TAX_FILE" > "${OUTDIR}/taxonomy_map.tsv" \
                              || cp -f "$TAX_FILE"  "${OUTDIR}/taxonomy_map.tsv"
  log "Taxonomy : $TAX_FILE ($(wc -l < "${OUTDIR}/taxonomy_map.tsv") genomes)"
fi

# =============================================================================
# STEP 1 — Extract nucleotide .fna files from archive
# Archive structure:
#   bac120_marker_genes_{reps|all}_r232/
#       {MARKER_ID}/
#           {genome_accession}_{MARKER_ID}.fna    (nucleotide)
#           {genome_accession}_{MARKER_ID}.faa    (amino acid)
# We extract only .fna (nucleotide) for KMA read mapping.
# =============================================================================
log "STEP 1: Extracting .fna files from archive"
EXTRACT_DIR="${OUTDIR}/.extracted"
mkdir -p "$EXTRACT_DIR"

t0=$(date +%s)
tar -xzf "$ARCHIVE" -C "$EXTRACT_DIR" \
    --wildcards '*.fna' --no-anchored 2>>"$LOG" \
  || die "tar extraction failed"
t1=$(date +%s)
n_fna=$(find "$EXTRACT_DIR" -name '*.fna' | wc -l)
log "  Extracted $n_fna .fna files in $(( t1-t0 ))s"

# Catalogue what we have
find "$EXTRACT_DIR" -type f -name '*.fna' -printf '%P\n' \
  | awk -F'/' '{fn=$NF; sub(/\.fna$/, "", fn); print fn}' \
  | sort -u > "${OUTDIR}/markers_present.txt"
n_markers=$(wc -l < "${OUTDIR}/markers_present.txt")
log "  Distinct markers : $n_markers (expect 120)"

find "$EXTRACT_DIR" -type f -name '*.fna' -print0 \
  | xargs -0 awk '/^>/ {sub(/^>/,""); print}' \
  | sort -u > "${OUTDIR}/genomes_present.txt"
n_genomes=$(wc -l < "${OUTDIR}/genomes_present.txt")
log "  Distinct genomes : $n_genomes"

# =============================================================================
# STEP 2 — Build combined FASTA with KMA/CCMetagen-compatible headers
#
# Header format:  >genome
#
# This format enables three levels of downstream analysis:
#   - Per-marker coverage  → extract file-basename before *.fna
#   - Per-genome coverage  → extract header from file
#   - Per-taxon abundance  → extract field after '|', split on ';'
#
# Example:
#   >RS_GCF_000001405.40
# =============================================================================
log "STEP 2: Concatenating FASTA with KMA-compatible headers"
RAW_FASTA="${OUTDIR}/gtdb_bac120_scg_${LABEL}.raw.fna"
HEADER_MAP="${OUTDIR}/header_map.tsv"

# Build awk taxonomy lookup table then process all .fna files
t0=$(date +%s)

find "$EXTRACT_DIR" -type f -name '*.fna' -print0 \
| python3 "$SCG_FASTA_BUILDER" \
    --taxonomy "${OUTDIR}/taxonomy_map.tsv" \
    --out-fasta "$RAW_FASTA" \
    --header-map "$HEADER_MAP"

t1=$(date +%s)
n_raw=$(grep -c '^>' "$RAW_FASTA" || echo 0)
log "  Raw sequences : $n_raw  ($(( t1-t0 ))s)"

# =============================================================================
# STEP 3 — Clean: length filter + deduplicate
# =============================================================================
CLEAN_FASTA="${OUTDIR}/gtdb_bac120_scg_${LABEL}.fna"

if [[ "$SKIP_CLEAN" == "1" ]]; then
  log "STEP 3: SKIP_CLEAN — using raw FASTA as final"
  cp -f "$RAW_FASTA" "$CLEAN_FASTA"
else
  log "STEP 3: Cleaning — length >= $MIN_LEN  +  deduplication"
  TMP_LEN="${OUTDIR}/.tmp_len.fna"
  TMP_DD="${OUTDIR}/.tmp_dedup.fna"

  seqkit seq  -j "$THREADS" -m "$MIN_LEN" "$RAW_FASTA" > "$TMP_LEN"
  n_len=$(grep -c '^>' "$TMP_LEN" || echo 0)
  log "  After length filter : $n_len"

  log "  Identifyng identical sequences: seqkit rmdup -s -i [exact sequence duplicates]"
  seqkit rmdup -j "$THREADS" -s -i "$TMP_LEN" \
    --id-regexp '^([^|]+\|[A-Za-z.]+(?:\s+[A-Za-z]+)?)' \
    -D "${OUTDIR}/duplicate_seqids.txt" \
    -d "${OUTDIR}/duplicates.txt" > "$TMP_DD"
  n_dd=$(grep -c '^>' "$TMP_DD" || echo 0)

  if [[ "$DEDUPLICATE_SEQUENCES" == "1" ]]; then
    cp -f "$TMP_DD" "$CLEAN_FASTA"
    log "  After dedup         : $n_dd  (removed $(( n_len - n_dd )) exact duplicates)"
  else
    cp -f "$TMP_LEN" "$CLEAN_FASTA"
  fi
  rm -f "$TMP_LEN" "$TMP_DD"
fi

pigz -p "$THREADS" -k -f "$CLEAN_FASTA"
n_final=$(grep -c '^>' "$CLEAN_FASTA" || echo 0)
log "  Final sequences  : $n_final"
log "  Plain FASTA      : $CLEAN_FASTA"
log "  Compressed       : ${CLEAN_FASTA}.gz ($(du -sh "${CLEAN_FASTA}.gz" | cut -f1))"

# =============================================================================
# STEP 4 — Marker × genome coverage report
# =============================================================================
log "STEP 4: Marker × genome coverage report"
COV="${OUTDIR}/marker_genome_coverage.tsv"
echo -e "marker\tn_sequences\tgenome_accessions" > "$COV"
grep '^>' "$CLEAN_FASTA" \
  | awk -F'[~|]' '
{
    # Remove leading ">" from genome ID
    genome = $1
    sub(/^>/, "", genome)
    # Extract marker ID
    marker = $2
    # Increment sequence count per marker
    count[marker]++
    # Append genome accession (semicolon-separated)
    if (marker in genomes) {
        genomes[marker] = genomes[marker] ";" genome
    } else {
        genomes[marker] = genome
    }
}
END {
    # Output: marker, count, genome list
    for (m in count)
        printf "%s\t%d\t%s\n", m, count[m], genomes[m]
    }
}
' \
| sort -k1,1 >> "$COV"
n_present=$(tail -n +2 "$COV" | wc -l)
log "  Markers in final FASTA: $n_present / 120"
[[ $n_present -lt 120 ]] && log "  [WARN] Missing markers — check $COV"

# =============================================================================
# STEP 5 — KMA indexing
# =============================================================================
KMA_DB="${OUTDIR}/gtdb_bac120_scg_${LABEL}"
if [[ "$HAS_KMA" == "1" ]]; then
  log "STEP 5: KMA indexing → $KMA_DB"
  t0=$(date +%s)
  kma index -i "$CLEAN_FASTA" -o "$KMA_DB" -verbose 2>&1 | tee -a "$LOG"
  t1=$(date +%s)
  log "  Index complete: $(( t1-t0 ))s"
  ls -lh "${KMA_DB}".* | tee -a "$LOG"
else
  log "STEP 5: KMA not available — run manually:"
  log "  kma index -i $CLEAN_FASTA -o $KMA_DB"
fi

# =============================================================================
# STEP 6 — Cleanup
# =============================================================================
if [[ "$KEEP_EXTRACTED" != "1" ]]; then
  log "STEP 6: Removing extracted intermediates"
  rm -rf "$EXTRACT_DIR" "$RAW_FASTA"
else
  log "STEP 6: KEEP_EXTRACTED=1 — retaining $EXTRACT_DIR"
fi

{
  echo "======================================================================"
  echo "COMPLETE: $(ts)"
  echo "KMA database prefix : $KMA_DB"
  echo "Header map          : $HEADER_MAP"
  echo "Marker coverage     : $COV"
  echo ""
  echo "Next — align reads:"
  echo "  kma -ipe R1.fq.gz R2.fq.gz -o sample_scg -t_db $KMA_DB \\"
  echo "      -1t1 -mem_mode -and -t $THREADS"
  echo ""
  echo "Next — compute GEQ in R:"
  echo "  Rscript gtdb_scg_geq.R --kma_res sample_scg.res \\"
  echo "      --header_map $HEADER_MAP --out sample_geq.tsv"
  echo "======================================================================"
} | tee -a "$LOG"

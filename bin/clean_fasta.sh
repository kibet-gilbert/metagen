#!/usr/bin/env bash
# =============================================================================
# clean_fasta.sh   —   MODULE 2 of 2
#
# Clean a raw FASTA file (from export_blastdb.sh OR download_all_refseq.sh)
# for use as a CCMetagen / KMA database:
#
#   Step 1 — Title filtering  (drop PREDICTED, vectors, synthetic, uncultured...)
#   Step 2 — Length filter    (keep sequences >= MIN_LEN)
#   Step 3 — Taxid-scoped dereplication  (exact + containment, per taxid)
#   Step 4 — Low-complexity masking (dustmasker → hard-mask to N)  [optional]
#   Step 5 — Compress + KMA index                                 [optional]
#
# STEP 3 (NEW) — taxid-scoped exact + containment dereplication
# -------------------------------------------------------------
# Sequence headers MUST begin with the numeric NCBI taxid as the first
# field, separated by HDR_DELIM (default '|'):   >taxID|description...
# Dereplication is performed INDEPENDENTLY within each taxid, so:
#   * identical sequences from DIFFERENT taxa are always retained;
#   * within a taxid, exact duplicates AND shorter sequences fully
#     contained (100% identity over 100% of their length) in a longer
#     one are collapsed, keeping the LONGEST as representative.
# This is done with mmseqs2 easy-linclust (run via a Singularity image),
# after partitioning the input by taxid. Singletons and any record whose
# header lacks a numeric taxid are passed through unchanged.
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
#   DEDUP_SEQUENCES=0        1=collapse dups/contained seqs; 0=report only (default)
#   COMPRESS_FINAL=1         1=gzip output (default); 0=plain
#   KEEP_INTERMEDIATE=0      1=keep tmp files; 0=delete (default)
#   THREADS=4                threads for seqkit / pigz / mmseqs
#
#   --- Step 3 (taxid-scoped dereplication) ---
#   HDR_DELIM='|'            header field separator; taxid is field 1
#   MMSEQS_SIF=<path.sif>    Singularity image providing `mmseqs` (REQUIRED
#                            when DEDUP_SEQUENCES=1)
#   SORT_MEM=8G              memory budget for the external sort (partitioning)
#   SORT_TMP=$OUTDIR/.sorttmp scratch dir for the external sort
# =============================================================================
set -Eeuo pipefail

on_error() {
  local exit_code=$?
  local line_no=$1
  echo "Error: exit code ${exit_code} at line ${line_no}" >&2
}

trap 'on_error $LINENO' ERR

# =============================================================================
# ── Defaults ──────────────────────────────────────────────────────────────────
# =============================================================================
MIN_LEN="${MIN_LEN:-300}"
MASK_LOW_COMPLEX="${MASK_LOW_COMPLEX:-0}"
DEDUP_SEQUENCES="${DEDUP_SEQUENCES:-0}"
COMPRESS_FINAL="${COMPRESS_FINAL:-1}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"
THREADS="${THREADS:-4}"

# Step 3 (taxid-scoped dereplication) settings
HDR_DELIM="${HDR_DELIM:-|}"
MMSEQS_SIF="${MMSEQS_SIF:-}"
SORT_MEM="${SORT_MEM:-8G}"

usage() {
  cat <<EOF
Usage: $0 <input.fasta[.gz]> [outdir]

  input   : raw FASTA (plain or .gz) — from export_blastdb.sh or download_all_refseq.sh
  outdir  : output directory (default: ./fasta_clean)

Environment:
  MIN_LEN=${MIN_LEN}
  MASK_LOW_COMPLEX=${MASK_LOW_COMPLEX}   (1=dustmasker, 0=skip)
  DEDUP_SEQUENCES=${DEDUP_SEQUENCES}     (1=collapse dups+contained, 0=report only)
  COMPRESS_FINAL=${COMPRESS_FINAL}       (1=gzip, 0=plain)
  KEEP_INTERMEDIATE=${KEEP_INTERMEDIATE} (1=keep tmp, 0=delete)
  THREADS=${THREADS}
  HDR_DELIM=${HDR_DELIM}                 (taxid is header field 1)
  MMSEQS_SIF=${MMSEQS_SIF:-<unset>}      (required when DEDUP_SEQUENCES=1)
EOF
}

INPUT="${1:-}"
OUTDIR="${2:-fasta_clean}"
[[ -z "$INPUT" ]] && { usage; exit 1; }
[[ ! -f "$INPUT" ]] && { echo "ERROR: Input not found: $INPUT"; exit 1; }
mkdir -p "$OUTDIR"
SORT_TMP="${SORT_TMP:-${OUTDIR}/.sorttmp}"

# =============================================================================
# ── Logging ───────────────────────────────────────────────────────────────────
# =============================================================================
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
HAS_KMA=1;    command -v kma    >/dev/null 2>&1 || { log "[WARN] kma absent — skipping index"; HAS_KMA=0; }
[[ "$MASK_LOW_COMPLEX" == "1" ]] && \
  command -v dustmasker >/dev/null 2>&1 || { log "[WARN] dustmasker not found — disabling masking"; MASK_LOW_COMPLEX=0; }
[[ "$COMPRESS_FINAL" == "1" ]] && \
  { PZIP=$(command -v pigz 2>/dev/null || command -v gzip 2>/dev/null || die "pigz/gzip not found"); }

# Step 3 needs a container runtime + an mmseqs image when actually collapsing.
SING=""
if [[ "$DEDUP_SEQUENCES" == "1" ]]; then
  SING=$(command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || true)
  [[ -n "$SING" ]] || die "DEDUP_SEQUENCES=1 but neither singularity nor apptainer found"
  [[ -n "$MMSEQS_SIF" && -f "$MMSEQS_SIF" ]] || die "DEDUP_SEQUENCES=1 requires MMSEQS_SIF=<path to mmseqs .sif>"
fi

# =============================================================================
# ── Helper: count sequences ───────────────────────────────────────────────────
# =============================================================================
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
KMA_DB="${OUTDIR}/${BASENAME}_ccmetagen"
[[ "$MASK_LOW_COMPLEX" == "1" ]] && FINAL="${OUTDIR}/${BASENAME}.ccmetagen.masked.fasta"
[[ "$COMPRESS_FINAL"   == "1" ]] && FINAL="${FINAL}.gz"

# ── Cleanup on exit ───────────────────────────────────────────────────────────
cleanup() {
  if [[ "$KEEP_INTERMEDIATE" != "1" ]]; then
    log "Cleaning intermediates in $TMP/"
    rm -rf "$TMP" "$SORT_TMP" || true
  else
    log "KEEP_INTERMEDIATE=1 — leaving $TMP/ intact"
  fi
}
trap cleanup EXIT

# =============================================================================
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

# =============================================================================
# ── STEP 1: Title filtering ───────────────────────────────────────────────────
# =============================================================================
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
  | grep -Eiv 'PREDICTED|vector|synthetic construct|cloning|patent|metagenome' \
  | awk 'BEGIN{RS=">"; ORS=""} NR>1 {print ">"$0}' \
  > "$STEP1"

n_STEP1=$(count_fa "$STEP1")
log "  After title filter: ${n_STEP1}"

# =============================================================================
# ── STEP 2: Length filter ─────────────────────────────────────────────────────
# =============================================================================
log "Step 2/4: Length filter (>= ${MIN_LEN} nt)"
#[WARN] you may switch on flag -g/--remove-gaps to remove spaces
seqkit seq -j "$THREADS" -m "$MIN_LEN" "$STEP1" > "$STEP2"
n_STEP2=$(count_fa "$STEP2")
log "  After length filter: ${n_STEP2}"

# =============================================================================
# ── STEP 3: Deduplicate ───────────────────────────────────────────────────────
# =============================================================================
# =============================================================================
# ── STEP 3: Taxid-scoped dereplication (exact + containment) ──────────────────
# =============================================================================
# Dereplicate WITHIN each taxid only. Headers must be >taxID<HDR_DELIM>...
# Approach:
#   3a. Partition: external-sort records by taxid, write one file per taxid.
#   3b. Per-taxid mmseqs easy-linclust (--min-seq-id 1.0 -c 1.0 --cov-mode 1):
#       collapses exact dups AND shorter contained seqs, keeping the longest.
#       Singletons / non-numeric-taxid records are passed through unchanged.
#   3c. Reassemble representatives + passthrough into STEP3.
# When DEDUP_SEQUENCES=0 we still PARTITION+REPORT how many would be removed,
# but pass STEP2 through unchanged (report-only, matching prior behaviour).
log "Step 3/4: Taxid-scoped dereplication (exact + containment, delim='${HDR_DELIM}')"

PARTS="${TMP}/parts"
REPS="${TMP}/reps"
mkdir -p "$PARTS" "$REPS" "$SORT_TMP"

# ---- 3a. Partition by taxid (bounded RAM via external sort) ------------------
# Linearise each record onto one line (newlines -> \x01), sort by taxid,
# then emit per-taxid files opening exactly one handle at a time.
COUNTS="${TMP}/taxid_counts.tsv"
: > "$COUNTS"
log "  3a. Partitioning by taxid (SORT_MEM=${SORT_MEM})"
awk -v delim="$HDR_DELIM" '
      BEGIN { RS=">"; ORS="" }
      NR==1 { next }
      {
        rec=$0; header=rec; sub(/\n.*/,"",header)
        split(header,a,delim); tax=a[1]
        if (tax !~ /^[0-9]+$/) tax="_NOTAX"
        enc=">" rec; gsub(/\n/,"\001",enc)
        print tax "\t" enc "\n"
      }
    ' "$STEP2" \
  | sort -t$'\t' -k1,1 -S "$SORT_MEM" -T "$SORT_TMP" --parallel="$THREADS" \
  | awk -v parts="$PARTS" -v countf="$COUNTS" '
      BEGIN { FS="\t"; cur=""; n=0 }
      {
        tax=$1; enc=$2
        if (tax!=cur) {
          if (cur!="") { print cur "\t" n >> countf; close(out) }
          cur=tax; n=0; out=parts "/" tax ".fasta"; printf "" > out
        }
        gsub(/\001/,"\n",enc); printf "%s", enc >> out; n++
      }
      END { if (cur!="") { print cur "\t" n >> countf; close(out) } }
    '
n_taxids=$(wc -l < "$COUNTS")
log "  3a. Partitioned into ${n_taxids} taxid groups"

# ---- 3b. Per-taxid dereplication --------------------------------------------
: > "$STEP3"
mm() { "$SING" exec --bind "$OUTDIR" "$MMSEQS_SIF" mmseqs "$@"; }

n_in_multi=0      # sequences entering multi-seq taxid dedup
n_rep_multi=0     # representatives kept from those
while IFS=$'\t' read -r tax cnt; do
  [[ -z "$tax" ]] && continue
  in="${PARTS}/${tax}.fasta"
  [[ -s "$in" ]] || continue

  # Passthrough: singletons and the non-numeric-taxid bucket are never collapsed
  if [[ "$tax" == "_NOTAX" || "$cnt" -le 1 ]]; then
    cat "$in" >> "$STEP3"
    continue
  fi

  if [[ "$DEDUP_SEQUENCES" == "1" ]]; then
    wt="${TMP}/mm_${tax}"; rm -rf "$wt"; mkdir -p "$wt"
    if mm easy-linclust "$in" "${wt}/clu" "${wt}/tmp" \
          --min-seq-id 1.0 -c 1.0 --cov-mode 1 \
          --threads "$THREADS" -v 1 >>"$LOG" 2>&1; then
      r=$(grep -c '^>' "${wt}/clu_rep_seq.fasta" 2>/dev/null || echo 0)
      cat "${wt}/clu_rep_seq.fasta" >> "$STEP3"
      n_in_multi=$(( n_in_multi + cnt ))
      n_rep_multi=$(( n_rep_multi + r ))
    else
      log "  [WARN] linclust failed for taxid ${tax} — passing through unchanged"
      cat "$in" >> "$STEP3"
      n_in_multi=$(( n_in_multi + cnt ))
      n_rep_multi=$(( n_rep_multi + cnt ))
    fi
    rm -rf "$wt"
  else
    # report-only: still measure collapse, but keep all sequences
    wt="${TMP}/mm_${tax}"; rm -rf "$wt"; mkdir -p "$wt"
    if mm easy-linclust "$in" "${wt}/clu" "${wt}/tmp" \
          --min-seq-id 1.0 -c 1.0 --cov-mode 1 \
          --threads "$THREADS" -v 1 >>"$LOG" 2>&1; then
      r=$(grep -c '^>' "${wt}/clu_rep_seq.fasta" 2>/dev/null || echo 0)
      n_in_multi=$(( n_in_multi + cnt ))
      n_rep_multi=$(( n_rep_multi + r ))
    fi
    rm -rf "$wt"
  fi
done < "$COUNTS"

n_STEP3=$(count_fa "$STEP3")
# Removable = (multi-seq inputs) - (representatives kept)
n_duplicates=$(( n_in_multi - n_rep_multi ))

if [[ "$DEDUP_SEQUENCES" == "1" ]]; then
  log "  Dereplicated within taxid (exact + containment)"
  log "  After dedup: ${n_STEP3} (removed ${n_duplicates} redundant sequences)"
  STEP4_INPUT=$STEP3
  n_STEP4_INPUT=${n_STEP3}
else
  log "  Dereplication skipped (reporting only)"
  log "  Would remove ${n_duplicates} redundant sequences across multi-seq taxids"
  STEP4_INPUT=$STEP2
  n_STEP4_INPUT=${n_STEP2}
fi

# =============================================================================
# ── STEP 4: Low-complexity masking ───────────────────────────────────────────
# =============================================================================
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
  n_STEP4=${n_STEP4_INPUT}
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

# =============================================================================
# STEP 5 — KMA indexing
# =============================================================================
if [[ "$HAS_KMA" == "1" ]]; then
  log "STEP 5: KMA indexing → $KMA_DB"
  t0=$(date +%s)
  kma index -i "$FINAL" -o "$KMA_DB" -verbose 2>&1 | tee -a "$LOG"
  t1=$(date +%s)
  log "  Index complete: $(( t1-t0 ))s"
  ls -lh "${KMA_DB}".* | tee -a "$LOG"
else
  log "STEP 5: KMA not available — run manually:"
  log ""
  log "Next step — build KMA index:"
  log "  kma index -i $FINAL -o $KMA_DB"
  log ""
  log "Done: $(ts)"
fi
# =============================================================================

{
  echo "======================================================================"
  echo "clean_fasta.sh SUMMARY"
  echo "Input sequences     : ${n_INPUT}"
  echo "After title filter  : ${n_STEP1}"
  echo "After length (${MIN_LEN}): ${n_STEP2}"
  echo "After dedup         : ${n_STEP4_INPUT}"
  echo "Redundant seqs      : ${n_duplicates}  (exact + contained, within-taxid)"
  echo "Final output        : $n_final sequences  ($sz_final)"
  echo "Output file         : $FINAL"
  echo "KMA database prefix : $KMA_DB"
  echo "Finished            : $(ts)"
  echo "======================================================================"
} | tee -a "$LOG"

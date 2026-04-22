#!/usr/bin/env bash
# =============================================================================
# download_all_refseq.sh
# Download NCBI RefSeq reference genomes, stamp taxID headers, audit,
# resolve unmapped genomes via NCBI eutils, merge to a single FASTA for
# CCMetagen/KMA, and optionally clean up intermediate files.
#
# Usage:
#   download_all_refseq.sh [workdir] [nthreads] [want_rep]
#
#   workdir  : destination folder          (default: refseq)
#   nthreads : parallel jobs               (default: 20)
#   want_rep : 1 = include representative  (default: 0 = reference only)
#
# Environment overrides:
#   KEEP_INTERMEDIATES=1   keep per-genome .fna.gz and .fna files after merge
#   SKIP_MERGE=1           skip the final cat/merge step
#   NCBI_EMAIL             e-mail for eutils polite requests (recommended)
# =============================================================================
set -Eeuo pipefail

# ── User settings ─────────────────────────────────────────────────────────────
workdir="${1:-refseq}"
threeds="${2:-20}"
want_rep="${3:-0}"
KEEP_INTERMEDIATES="${KEEP_INTERMEDIATES:-0}"
SKIP_MERGE="${SKIP_MERGE:-0}"
NCBI_EMAIL="${NCBI_EMAIL:-${USER}@ilri.org}"

[[ -z "$workdir" ]] && { echo "[ERROR] Provide a destination folder"; exit 1; }
[[ -z "$threeds" ]] && threeds=20
mkdir -p "$workdir"

# ── Logging ───────────────────────────────────────────────────────────────────
LOG="${workdir}/download_refseq.log"
ts()  { date +"%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*" | tee -a "$LOG"; }
die() { echo "[$(ts)] ERROR: $*" | tee -a "$LOG" >&2; exit 1; }

{
  echo "======================================================================"
  echo "Started          : $(ts)"
  echo "Script           : $0"
  echo "workdir          : $workdir"
  echo "Threads          : $threeds"
  echo "want_rep         : $want_rep"
  echo "KEEP_INTERMEDIATES: $KEEP_INTERMEDIATES"
  echo "SKIP_MERGE       : $SKIP_MERGE"
  echo "======================================================================"
} | tee -a "$LOG"

# ── Tool checks ───────────────────────────────────────────────────────────────
for t in wget parallel pigz gzip awk sed python3 curl; do
  command -v "$t" >/dev/null 2>&1 || die "Missing required tool: $t"
done

# =============================================================================
# STEP 1 — Download assembly summary
# =============================================================================
log "STEP 1: Download NCBI RefSeq assembly summary"
refseq_assembly_summary="${workdir}/assembly_summary.txt"
if [[ -s "$refseq_assembly_summary" ]]; then
  log "  Found existing: $refseq_assembly_summary — skipping download"
else
  wget -q -O "$refseq_assembly_summary" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt" \
    2>>"$LOG" || die "Failed to download assembly summary"
  log "  Downloaded: $refseq_assembly_summary"
fi
log "  Entries: $(wc -l < "$refseq_assembly_summary")"

# =============================================================================
# STEP 2 — Filter to reference (and optionally representative) genome URLs
# =============================================================================
log "STEP 2: Filter assembly summary → FTP URL list"
assembly_summary_ref="${workdir}/assembly_summary_ref.txt"
if [[ -s "$assembly_summary_ref" ]]; then
  log "  Found existing: $assembly_summary_ref — skipping generation"
else
  if [[ "$want_rep" == "1" ]]; then
    awk -F '\t' 'NR>1 && $1 !~ /^#/ && ($5=="reference genome" || $5=="representative genome") {print $20}' \
      "$refseq_assembly_summary" \
    | sed -e 's/^"//; s/"$//' -e 's:/*$::' \
    > "$assembly_summary_ref"
  else
    awk -F '\t' 'NR>1 && $1 !~ /^#/ && $5=="reference genome" {print $20}' \
      "$refseq_assembly_summary" \
    | sed -e 's/^"//; s/"$//' -e 's:/*$::' \
    > "$assembly_summary_ref"
  fi
  log "  Generated: $assembly_summary_ref"
fi
refseq_URLs="$assembly_summary_ref"
log "  URLs to download: $(wc -l < "$refseq_URLs")"

# =============================================================================
# STEP 3 — Parallel download with retries
# =============================================================================
log "STEP 3: Download genome FASTA files in parallel"

_download_job() {
  set -Eeuo pipefail
  local idx="$1"
  local url="${2%/}"
  local base
  base="$(basename "$url")"
  local file="$WORKDIR/${base}_genomic.fna.gz"
  local full="$url/${base}_genomic.fna.gz"
  local md5url="$url/md5checksums.txt"
  local md5check_file="${WORKDIR}/_dl/md5checksum.r${ROUND}.txt"

  if [[ -s "$file" ]]; then
    echo "STATUS:OK $url" >&2
  else
    echo "STATUS:GET $full" >&2
    local tmp="${file}.part.$$"
    if wget -q \
            --user-agent="ILRI-ccmetagen/1.0 (+contact: ${NCBI_EMAIL:-user@example.org})" \
            --tries=3 --timeout=90 --read-timeout=600 \
            --waitretry=5 --retry-on-http-error=429,500,502,503,504 \
            -O "$tmp" "$full"; then
      mv "$tmp" "$file"
      echo "STATUS:DONE $url" >&2
    else
      rm -f "$tmp" || true
      echo "STATUS:FAIL $url" >&2
      return 1
    fi
  fi

  # MD5 capture (skip if already recorded this round)
  local md5sumcheck md5raw md5line
  md5sumcheck=$(grep "${base}_genomic.fna.gz" "$md5check_file" 2>/dev/null || true)
  if [[ -n "$md5sumcheck" ]]; then return 0; fi

  md5raw="$(wget -q -O - "$md5url" 2>/dev/null || true)"
  if [[ -n "$md5raw" ]]; then
    md5line="$(printf '%s\n' "$md5raw" | grep "${base}_genomic.fna.gz" || true)"
    if [[ -n "$md5line" ]]; then
      printf '%s\n' "$md5line" >> "$md5check_file"
    fi
  fi
}
export -f _download_job
export NCBI_EMAIL

download_refseq_polite() {
  local workdir="$1"  max_rounds="$2"  jobs="$3"  delay="$4"  src="$5"
  local state="${workdir}/_dl"
  mkdir -p "$state"
  export WORKDIR="$workdir"

  sed -E -i \
    -e "s/^[\"']//" -e "s/[\"']$//" -e 's:/*$::' \
    "$src"

  cp -f "$src" "${state}/urls.r1.txt"

  local r=1  next_fail=""
  while (( r <= max_rounds )); do
    local in="${state}/urls.r${r}.txt"
    if [[ ! -s "$in" ]]; then
      log "  Nothing to do in round ${r} — all done."
      break
    fi

    local log_r="${state}/wget.r${r}.log"
    local md5f="${state}/md5checksum.r${r}.txt"
    local jlog="${state}/parallel.r${r}.joblog"
    : > "$log_r"; : > "$md5f"

    local total
    total=$(wc -l < "$in")
    log "  Round ${r}/${max_rounds} — ${total} URLs"

    nl -ba -w1 -s $'\t' "$in" > "${state}/urls.r${r}.num.txt"

    export ROUND="$r"
    export MAX_ROUNDS="$max_rounds"
    export TOTAL="$total"
    export LOGFILE="$log_r"

    parallel --jobs "$jobs" --delay "$delay" \
      --colsep '\t' --silent \
      --joblog "$jlog" \
      _download_job {1} {2} \
      :::: "${state}/urls.r${r}.num.txt" >> "$log_r" 2>&1 || true

    next_fail="${state}/urls.r$((r+1)).txt"
    awk '/^STATUS:FAIL /{print $2}' "$log_r" | sort -u > "$next_fail"

    local n_get n_done n_ok n_fail
    n_get=$(grep -c '^STATUS:GET'  "$log_r" || true)
    n_done=$(grep -c '^STATUS:DONE' "$log_r" || true)
    n_ok=$(grep -c  '^STATUS:OK'   "$log_r" || true)
    n_fail=$(grep -c '^STATUS:FAIL' "$log_r" || true)
    log "  Round ${r} summary: GET=${n_get} DONE=${n_done} OK=${n_ok} FAIL=${n_fail}"

    if [[ ! -s "$next_fail" ]]; then
      log "  All downloads complete after round ${r}."
      break
    fi

    local sleep_s=$(( 5 * (2 ** (r-1)) ))
    (( sleep_s > 300 )) && sleep_s=300
    log "  Sleeping ${sleep_s}s before round $((r+1)) (rate-limit backoff)..."
    sleep "$sleep_s"
    (( r++ ))
  done

  if [[ -n "$next_fail" && -s "$next_fail" ]]; then
    cp -f "$next_fail" "${workdir}/download_failures.final.txt"
    log "  [WARN] Persistent failures: ${workdir}/download_failures.final.txt"
  else
    log "  No remaining failures."
  fi
}

download_refseq_polite "$workdir" 10 "$threeds" 0.25 "$refseq_URLs"

log "  Downloaded .fna.gz count: $(ls "$workdir"/*_genomic.fna.gz 2>/dev/null | wc -l)"

# =============================================================================
# STEP 4 — Download accession-to-taxid maps
# =============================================================================
log "STEP 4: Download accession-to-taxid maps"

for pair in \
    "nucl_gb.accession2taxid.gz|https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz" \
    "nucl_wgs.accession2taxid.gz|https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz"; do
  fname="${pair%%|*}"
  url="${pair##*|}"
  dest="${workdir}/${fname}"
  if [[ -s "$dest" ]]; then
    log "  Found existing: $dest — skipping"
  else
    log "  Downloading: $dest"
    wget -q -O "$dest" "$url" 2>>"$LOG" || log "  [WARN] Failed to download $url"
  fi
done

log "STEP 4b: Build versioned_accession → taxid map"
acc_taxid_map="${workdir}/accession_taxid_nucl-rs.map"
if [[ -s "$acc_taxid_map" ]]; then
  log "  Found existing: $acc_taxid_map — skipping"
else
  zgrep "_" "$workdir"/nucl_*.accession2taxid.gz \
    | cut -f2-3 > "$acc_taxid_map"
  log "  Map entries: $(wc -l < "$acc_taxid_map")"
fi

# =============================================================================
# STEP 5 — Build accession → filepath map
# =============================================================================
log "STEP 5: Build accession → filepath map"
acc_fqpath_map="${workdir}/accession_fqpath-rs.map"
if [[ -s "$acc_fqpath_map" ]]; then
  log "  Found existing: $acc_fqpath_map — skipping"
else
  accession_filename() {
    zcat "$1" | awk -v f="$1" '/^>/{print substr($1,2) "\t" f}'
  }
  export -f accession_filename
  parallel --keep-order -j "$threeds" \
    accession_filename {} \
    ::: "$workdir"/*_genomic.fna.gz \
    > "$acc_fqpath_map"
  log "  Accession-path pairs: $(wc -l < "$acc_fqpath_map")"
fi

# =============================================================================
# STEP 6 — Join maps → taxid_fqpath.map
# =============================================================================
log "STEP 6: Join accession maps → taxid_fqpath.map"
taxid_fqpath="${workdir}/taxid_fqpath.map"
if [[ -s "$taxid_fqpath" ]]; then
  log "  Found existing: $taxid_fqpath — skipping"
else
  join -j1 -e MISSING -o 1.2,2.2 \
    <(sort "$acc_taxid_map") \
    <(sort "$acc_fqpath_map") \
    | tr ' ' '\t' \
    | sed 's|\t\./|\t|' \
    | sed 's|//|/|g' \
    | sort -u > "$taxid_fqpath"
  log "  taxid-filepath entries: $(wc -l < "$taxid_fqpath")"
fi

# =============================================================================
# STEP 7 — Stamp taxID headers (parallel, single-pass via pigz)
# =============================================================================
log "STEP 7: Stamp taxID into FASTA headers"

_stamp_job() {
  local taxid="$1"
  local src="$2"
  local base
  base="$(basename "${src%_genomic.fna.gz}")"
  local out="${src%_genomic.fna.gz}_genomic_taxID2.fna.gz"

  if [[ -s "$out" ]]; then
    echo "[INFO] Skipping ${base} — output exists"
    return 0
  fi

  gzip -dc "$src" \
    | sed -E "s/^>/>$taxid|/" \
    | pigz -p 2 -c \
    > "$out" \
  && echo "[INFO] Stamped: $(basename "$out")" \
  || { echo "[ERROR] Failed: $src" >&2; return 1; }
}
export -f _stamp_job

parallel --colsep '\t' -j "$threeds" \
  _stamp_job {1} {2} \
  :::: "$taxid_fqpath" \
  >> "${workdir}/stamplog.temp" 2>&1

log "  Stamped .fna.gz count: $(ls "$workdir"/*_taxID2.fna.gz 2>/dev/null | wc -l)"

# =============================================================================
# STEP 8 — Audit: cross-reference .fna.gz, taxid_fqpath.map, and stamped outputs
# =============================================================================
log "STEP 8: Audit — cross-reference files vs maps"
audit_file="${workdir}/audit.txt"

awk -F'\t' '
  NR==FNR {
    n = split($1, a, "/"); on_disk[a[n]] = $1; next
  }
  {
    path = $2; gsub(/\/\//, "/", path)
    n = split(path, a, "/"); base = a[n]
    in_map[base] = path; taxid[base] = $1
  }
  END {
    for (b in in_map) if (!(b in on_disk)) print "MAP_NO_FILE\t" taxid[b] "\t" in_map[b]
    for (b in on_disk) if (!(b in in_map)) print "FILE_NO_MAP\t\t" on_disk[b]
  }
' <(ls "$workdir"/*_genomic.fna.gz 2>/dev/null) "$taxid_fqpath" \
  | sort > "$audit_file"

n_map_no_file=$(grep -c '^MAP_NO_FILE' "$audit_file" || true)
n_file_no_map=$(grep -c '^FILE_NO_MAP' "$audit_file" || true)
log "  In map but missing on disk  : $n_map_no_file"
log "  On disk but missing from map: $n_file_no_map"
cat "$audit_file" >> "$LOG"

# =============================================================================
# STEP 9 — Resolve unmapped genomes via NCBI eutils
# =============================================================================
log "STEP 9: Resolve unmapped genomes via NCBI eutils"

missing_gcf="${workdir}/missing_gcf_accessions.txt"
missing_taxids="${workdir}/missing_taxids.txt"
residual_map="${workdir}/taxid_fqpath.residual_fs.map"
: > "$residual_map"

if [[ "$n_file_no_map" -gt 0 ]]; then
  # Extract GCF accession from FILE_NO_MAP lines
  awk '$1=="FILE_NO_MAP"{print $3}' "$audit_file" \
    | sed 's|.*/||; s/_[^_]*_genomic\.fna\.gz$//' \
    | sort -u > "$missing_gcf"
  log "  Unmapped GCF accessions: $(wc -l < "$missing_gcf")"

  # Peek at accessions inside unmapped files
  log "  Checking accessions inside unmapped files:"
  while IFS= read -r f; do
    echo "=== $(basename "$f") ===" >> "$LOG"
    zcat "$f" | grep '^>' | head -3 >> "$LOG"
  done < <(awk '$1=="FILE_NO_MAP"{print $3}' "$audit_file")

  # Fetch taxids from NCBI eutils
  log "  Fetching taxids from NCBI eutils (polite: 0.4s delay)..."
  : > "$missing_taxids"
  while IFS= read -r gcf; do
    uid=$(curl -s \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=${gcf}&retmode=json&email=${NCBI_EMAIL}" \
      | python3 -c "
import sys,json
d=json.load(sys.stdin)
ids=d.get('esearchresult',{}).get('idlist',[])
print(ids[0] if ids else '')
" 2>/dev/null)

    if [[ -z "$uid" ]]; then
      echo -e "NA\t$gcf" >> "$missing_taxids"
      sleep 0.4; continue
    fi

    tax=$(curl -s \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=${uid}&retmode=json&email=${NCBI_EMAIL}" \
      | python3 -c "
import sys,json
d=json.load(sys.stdin)
r=d.get('result',{})
u=r.get('uids',[])
print(r[u[0]].get('taxid','NA') if u else 'NA')
" 2>/dev/null)

    echo -e "$tax\t$gcf" >> "$missing_taxids"
    sleep 0.4
  done < "$missing_gcf"

  log "  Taxids resolved: $(grep -v '^NA' "$missing_taxids" | wc -l) / $(wc -l < "$missing_taxids")"

  # Build residual taxid → filepath entries
  while IFS=$'\t' read -r taxid gcf; do
    [[ "$taxid" == "NA" || -z "$taxid" ]] && {
      log "  [WARN] No taxid for $gcf — skipping"
      continue
    }
    fpath=$(ls "$workdir/${gcf}_"*_genomic.fna.gz 2>/dev/null | head -1 || true)
    if [[ -n "$fpath" ]]; then
      echo -e "$taxid\t$fpath" >> "$residual_map"
    else
      log "  [WARN] File not found for $gcf"
    fi
  done < "$missing_taxids"

  log "  Residual map entries: $(wc -l < "$residual_map")"

  # Stamp the residual genomes
  if [[ -s "$residual_map" ]]; then
    log "  Stamping residual genomes..."
    parallel --colsep '\t' -j "$threeds" \
      _stamp_job {1} {2} \
      :::: "$residual_map" \
      >> "${workdir}/stamplog.residual.temp" 2>&1
  fi
else
  log "  No unmapped genomes — skipping eutils lookup."
fi

# =============================================================================
# STEP 10 — Final stamp audit
# =============================================================================
log "STEP 10: Final stamp audit"
n_fna_gz=$(ls "$workdir"/*_genomic.fna.gz 2>/dev/null | wc -l || echo 0)
n_stamped=$(ls "$workdir"/*_taxID2.fna.gz  2>/dev/null | wc -l || echo 0)
n_unstamped=$(( n_fna_gz - n_stamped ))
log "  Input  .fna.gz    : $n_fna_gz"
log "  Stamped _taxID2.fna.gz: $n_stamped"
log "  Unstamped          : $n_unstamped"

# Find any remaining unstamped files
unstamped_list="${workdir}/unstamped.txt"
for f in "$workdir"/*_genomic.fna.gz; do
  expected="${f%_genomic.fna.gz}_genomic_taxID2.fna.gz"
  [[ ! -s "$expected" ]] && echo "$f"
done > "$unstamped_list"
n_truly_missing=$(wc -l < "$unstamped_list")
log "  Still missing stamped output: $n_truly_missing"
if [[ $n_truly_missing -gt 0 ]]; then
  log "  [WARN] See: $unstamped_list"
fi

# =============================================================================
# STEP 11 — Merge all stamped FASTA into one file for CCMetagen/KMA
# =============================================================================
if [[ "$SKIP_MERGE" == "1" ]]; then
  log "STEP 11: SKIP_MERGE=1 — skipping merge"
else
  log "STEP 11: Merge all *_taxID2.fna.gz into single FASTA for KMA indexing"
  merged="${workdir}/refseq_ccmetagen.fna.gz"

  if [[ -s "$merged" ]]; then
    log "  Found existing merged file: $merged — skipping"
  else
    n_to_merge=$(ls "$workdir"/*_taxID2.fna.gz 2>/dev/null | wc -l || echo 0)
    log "  Merging $n_to_merge files → $merged"
    t0=$(date +%s)
    cat "$workdir"/*_taxID2.fna.gz > "$merged"
    t1=$(date +%s)
    merged_size=$(du -sh "$merged" | cut -f1)
    log "  Merge complete: ${merged_size} in $(( t1-t0 ))s"
    log "  Sequence count: $(zgrep -c '^>' "$merged")"
  fi
fi

# =============================================================================
# STEP 12 — Cleanup intermediates
# =============================================================================
if [[ "$KEEP_INTERMEDIATES" == "1" ]]; then
  log "STEP 12: KEEP_INTERMEDIATES=1 — retaining all files"
else
  log "STEP 12: Cleaning intermediate files"

  # Verify merge exists and is non-empty before removing sources
  merged="${workdir}/refseq_ccmetagen.fna.gz"
  if [[ ! -s "$merged" && "$SKIP_MERGE" != "1" ]]; then
    log "  [WARN] Merged file missing or empty — NOT removing intermediates"
  else
    # Remove per-genome source .fna.gz (originals)
    log "  Removing *_genomic.fna.gz source files..."
    find "$workdir" -maxdepth 1 -name '*_genomic.fna.gz' \
      ! -name '*_taxID2*' -delete
    log "  Removing *_taxID2.fna.gz (already merged)..."
    find "$workdir" -maxdepth 1 -name '*_taxID2.fna.gz' -delete
    # Remove per-round download logs (keep summary log)
    log "  Removing per-round download logs..."
    rm -rf "${workdir}/_dl" || true
    log "  Cleanup complete. Retained:"
    log "    $merged        (merged CCMetagen FASTA)"
    log "    $taxid_fqpath  (taxid map)"
    log "    $acc_taxid_map (accession map)"
    log "    $LOG           (this log)"
    log "    ${workdir}/audit.txt"
    log "    ${workdir}/download_failures.final.txt (if any)"
    log "    ${workdir}/missing_*.txt (if any)"
  fi
fi

# =============================================================================
log "======================================================================"
log "DONE: $(ts)"
log "  Merged FASTA for KMA: ${workdir}/refseq_ccmetagen.fna.gz"
log "  Full log: $LOG"
log "======================================================================"

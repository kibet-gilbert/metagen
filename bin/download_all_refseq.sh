#!/usr/bin/env bash

# This script was updated from: https://github.com/vrmarcelino/CCMetagen/blob/master/docs/benchmarking/rename_nt/rename_refseq.sh
#=========================================================================#
# download_refseq.sh does the following:
# 1. Download all "reference genomes" from NCBI Refseq, accessions and taxids
# 2. Edit headers to >taxid|accession.
# NOTE: For CCMetagen.py add the '-r nt' command, ('-r nr' is the default)

# Usage: download_all_refseq.sh  [$dest_folder]  [$nthreads] [$want_rep]
#--------------------------------------------------------------------------#

#!/usr/bin/env bash
set -Eeuo pipefail

# ---- user settings ----
# Assumptions: threeds=20, workdir="refseq".
workdir="${1:-refseq}"
threeds="${2:-20}"     # number of parallel downloads
want_rep="${3:-0}"    # set to 1 to include "representative genome" in addition to "reference genome"
# -----------------------

if [ -z "$workdir" ];
then echo " < ! >   need to provide a destination folder!" ;
exit 1 ; fi

if [ -z "$threeds" ] ; then threeds=10 ; fi

# the dir wherin the RefSeq is stowed
mkdir -p $workdir
echo -e "\n[Workdir]: ${workdir}"
echo -e "[Threads]: ${threeds}"

echo " +-1-+ recommended to download NCBI RefSeq genomes - 30m  ---------"
# 1) Download the latest assembly summary for RefSeq
refseq_assembly_summary="${workdir}/assembly_summary.txt"
if [[ -s "$refseq_assembly_summary" ]]; then
    echo "[INFO] Found existing assembly summary, skipping download: $refseq_assembly_summary"
else
    echo "[INFO] Downloading NCBI RefSeq assembly summary …"
    wget -O "$refseq_assembly_summary" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt" \
        > "$workdir/wgetlog1.temp" 2>&1
fi
echo -e "\t[NCBI RefSeq assembly summary file]: ${refseq_assembly_summary}"

# 2) Filter lines - clean URL list:
#    - keep refseq_category == "reference genome" (and optionally "representative genome")
#    - extract column 20 (ftp_path)
# Column 20 is ftp_path; column 5 is refseq_category.
assembly_summary_ref="${workdir}/assembly_summary_ref.txt"
if [[ -s "$assembly_summary_ref" ]]; then
    echo "[INFO] Found existing assembly URLs file, skipping generation: $assembly_summary_ref"
else
    echo "[INFO] Generating NCBI RefSeq assembly summary URLs…"
    if [[ "${want_rep}" == "1" ]]; then
      awk -F '\t' 'NR>1 && $1 !~ /^#/ && ($5=="reference genome" || $5=="representative genome") {print $20}' \
        "${refseq_assembly_summary}" \
      | sed -e 's/^"//; s/"$//' -e 's:/*$::' \
      > "${assembly_summary_ref}"
    else
      awk -F '\t' 'NR>1 && $1 !~ /^#/ && $5=="reference genome" {print $20}' \
        "${refseq_assembly_summary}" \
      | sed -e 's/^"//; s/"$//' -e 's:/*$::' \
      > "${assembly_summary_ref}"
    fi
fi

refseq_URLs="${assembly_summary_ref}"
echo -e "\t[NCBI RefSeq genome assembly urls]: ${refseq_URLs}"

# 3) Download in parallel, robustly
#    - construct basename safely after trimming trailing slash
#    - retries + resume (-c)
#    - write to ${workdir}/${basename}_genomic.fna.gz
#    - quiet progress but keep errors in log
# Polite RefSeq downloader with retries
# Usage: download_refseq_polite refseq 10 8 0.25
#   arg1 = workdir (defaults: refseq)
#   arg2 = max_rounds (defaults: 10)
#   arg3 = jobs (defaults: 8)
#   arg4 = delay seconds between job starts (defaults: 0.25)
#   arg5 = NCBI RefSeq genome assembly urls (defaults: ${refseq_URLs})
# Define the per-job worker as a real function and export it
_download_job() {
  set -Eeuo pipefail
  local idx="$1"
  local url="${2%/}"
  local base="$(basename "$url")"
  local file="$WORKDIR/${base}_genomic.fna.gz"
  local full="$url/${base}_genomic.fna.gz"
  local md5url="$url/md5checksums.txt"
  local md5check_file="${WORKDIR}/_dl/md5checksum.r${ROUND}.txt"

  # Skip if $file is already present and non-empty
  if [[ -s "$file" ]]; then
    echo "[INFO] Round ${ROUND}/${MAX_ROUNDS}: ${base} — ${idx}/${TOTAL} skipping..." #| tee -a "$LOGFILE"
    echo "STATUS:OK $url" >&2
    # return 0
  else
    echo "[INFO] Round ${ROUND}/${MAX_ROUNDS}: ${base} — ${idx}/${TOTAL} fetching..." #| tee -a "$LOGFILE"
    echo "STATUS:GET $full" >&2
    local tmp="${file}.part.$$"

    if wget -q \
            --user-agent="ILRI-ccmetagen/1.0 (+contact: ${USER}@ilri.org)" \
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

  # Only reached after a successful fresh download 
  # Check if md5 was already recorded from a previous attempt this round
  # MD5 capture
  local md5raw md5line md5sumcheck
  md5sumcheck=$(grep "${base}_genomic.fna.gz" "${md5check_file}" || true)
  if [[ -n "$md5sumcheck" ]]; then
    return 0
  fi

  echo "[INFO] md5sumcheck: ${base} — ${idx}/${TOTAL} fetching..."
  md5raw="$(wget -q -O - "$md5url" 2>/dev/null || true)"
  if [[ -n "$md5raw" ]]; then
    md5line="$(printf '%s\n' "$md5raw" | grep "${base}_genomic.fna.gz" || true)"
    if [[ -n "$md5line" ]]; then
      printf '%s\t%s\n' "$md5line" >> "${md5check_file}"
    else
      echo -e "\t[INFO]: NO_MD5_LINE $url"
    fi
  else
    echo -e "\t[INFO]: NO_MD5_FILE $url"
  fi
}
export -f _download_job

download_refseq_polite() {
  set -Eeuo pipefail
  local workdir="${1:-refseq}"
  local max_rounds="${2:-10}"
  local jobs="${3:-${threeds}}"
  local delay="${4:-0.25}"
  local src="${5:-${refseq_URLs}}"
  local state="${workdir}/_dl"
  mkdir -p "$state"

  # Ensure URL list is clean (no quotes, no trailing slashes)
  sed -E -i \
    -e "s/^[\"']//" \
    -e "s/[\"']$//" \
    -e 's:/*$::' \
    "$src"

  # Seed input list (round 1 uses the full list)
  cp -f "$src" "${state}/urls.r1.txt"

  local r=1
  while (( r <= max_rounds )); do
    local in="${state}/urls.r${r}.txt"
    # If there's nothing to download, stop early
    if [[ ! -s "$in" ]]; then
      echo "[INFO] Nothing to do in round ${r}—all done."
      break
    fi

    # If there is - set per‑round paths & clears logs
    local log="${state}/wget.r${r}.log"
    local md5CheckFile="${state}/md5checksum.r${r}.txt"
    local joblog="${state}/parallel.r${r}.joblog"
    : > "$log"
    : > "$md5CheckFile"
    local total=$(wc -l < "$in")
    echo "[INFO] wget logs: ${log}..."
    echo "[INFO] Parallel logs: ${joblog}..."
    echo "[INFO] md5checksum file: ${md5CheckFile}..."
    echo "[INFO] Round ${r}/${max_rounds} — ${total} URLs to fetch..."

    # Build a numbered input for this round: <idx>\t<url>
    nl -ba -w1 -s $'\t' "$in" > "${state}/urls.r${r}.num.txt"
    local urlNumFile="${state}/urls.r${r}.num.txt"
    
    # Export metadata so jobs can log “Round R/M : BASE — idx/total <status>”
    export WORKDIR="$workdir"
    export LOGFILE="$log"
    export ROUND="$r"
    export MAX_ROUNDS="$max_rounds"
    export TOTAL="$total"
    local next_fail=""

    # Avoid --halt now,fail=1 here—we want all jobs to run, then we retry the remaining failures as a batch.
    parallel --jobs "$jobs" --delay "$delay" \
         --colsep '\t' --silent \
         --joblog "$joblog" \
         _download_job {1} {2} \
         :::: "$urlNumFile" >> "$log" 2>&1 || true

    # Build next-round retry list from FAIL lines
    # save the failure path after building it:
    local next_fail="${state}/urls.r$((r+1)).txt"
    awk '/^FAIL /{print $2}' "$log" | sort -u > "$next_fail"

    # Progress summary
    printf "[INFO] Round %d summary: " "$r"
    awk '
      BEGIN{ok=0;done=0;fail=0;get=0}
      /^STATUS:OK/ {ok++}
      /^STATUS:DONE/ {done++}
      /^STATUS:FAIL/ {fail++}
      /^STATUS:GET/ {get++}
      END{
        printf("GET=%d DONE=%d OK=%d FAIL=%d\n", get, done, ok, fail)
      }' "$log"

    # Stop if nothing failed
    if [[ ! -s "$next_fail" ]]; then
      echo "[INFO] All downloads complete after round ${r}."
      break
    fi

    # Optional: exponential backoff between rounds (5s, 10s, 20s, ...)
    local sleep_s=$(( 5 * (2 ** (r-1)) ))
    # Cap the sleep to something reasonable (e.g. 300s) in long runs
    (( sleep_s > 300 )) && sleep_s=300
    echo "[INFO] Sleeping ${sleep_s}s before next round \(to ease 503/429 throttling\) ..."
    sleep "$sleep_s"
    (( r++ ))
  done

  # Write a final failure list (if any)
  if [[ -s "$next_fail" ]]; then
    cp -f "$next_fail" "${workdir}/download_failures.final.txt"
    echo "[WARN] Some URLs still failed after ${max_rounds} rounds:"
    echo "       see ${workdir}/download_failures.final.txt"
  else
    echo "[INFO] No remaining failures."
  fi

  echo "[INFO] Logs: ${state}/wget.r*.log; joblogs: ${state}/parallel.r*.joblog"
}

# Run with defaults (workdir=refseq, rounds=10, jobs=8, delay=0.25s)
download_refseq_polite "$workdir" 10 "$threeds" 0.25 "${refseq_URLs}"

echo " + 2 +   get RefSeq ( _ ) accession-taxid listings - (large) 10m?------"
# These files are very large.
# accessions you expect in FASTA deflines are in nucl_gb and nucl_wgs
wget -O "$workdir/nucl_gb.accession2taxid.gz" \
  https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz >> $workdir/wgetlog2.temp 2>&1
# nucl_wgs
wget -O "$workdir/nucl_wgs.accession2taxid.gz" \
  https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz >> $workdir/wgetlog3.temp 2>&1
# Build a single versioned_accession -> taxid map from the gz files
zgrep "_" $workdir/nucl_*.accession2taxid.gz \
  | cut -f 2-3 > $workdir/accession_taxid_nucl-rs.map

# Build "accession" → "file" map from your downloaded FASTAs
echo " + 3 +    list accession -> fasta-file pairs - (stream) 3m  --------"
# We want: accession<TAB>file_path
accession_filename () {
  # $1 = gz fasta file path
  zcat $1 | awk -v filen=$1 '/^>/ {print substr($1, 2), "\t", filen}'
}
export -f accession_filename

# Run in parallel over all *genomic.fna.gz
# (Keep order optional; we just want all pairs.)
# threeds=20
parallel --keep-order -j $threeds \
   accession_filename {} \
   ::: "${workdir}"/*genomic.fna.gz > "$workdir/accession_fqpath-rs.map"

# Join maps (versioned accession) → taxid + path
echo " + 4 +   sort & join two lists by accession - 3m ---(ignore 9606 error)-"
# Left join on accession.version (field 1 in both maps)
# We want: taxid<TAB>file_path
time join -j 1 -e XXXcatastrophotronXXXX -o 1.2,2.2 \
  <( sort "$workdir/accession_taxid_nucl-rs.map" ) \
  <( sort "$workdir/accession_fqpath-rs.map" ) \
  | tr " " "\t" \
  | sort -u > "$workdir/taxid_fqpath.map"

# Re‑stamp headers to >taxid|description (idempotent, streaming)
echo " + 5 +   re-stamp headers - 30m  --------------------------------"

_stamp_job() {
  local taxid="$1"
  local src="$2"
  local base
  base="$(basename "${src%.fna.gz}")"
  local out="${src%_genomic.fna.gz}_genomic_taxID2.fna.gz"

  # Skip if already present and non-empty
  if [[ -s "$out" ]]; then
    echo "[INFO] Skipping ${base} — output already exists."
    return 0
  fi

  echo "[INFO] Stamping ${base} with taxid ${taxid}..."
  gzip -dc "$src" \
    | sed -E "s/^>/>$taxid|/" \
    | pigz -p 2 -c \
    > "$out" \
  && echo "[INFO] Done: $(basename "$out")" \
  || { echo "[ERROR] Failed: $src" >&2; return 1; }
}
export -f _stamp_job

parallel --colsep "\t" -j "$threeds" \
  _stamp_job {1} {2} \
  :::: "$workdir/taxid_fqpath.map" \
  >> "$workdir/stamplog.temp" 2>&1

echo " +-6-+   done. NOTE :: use the '-r nt' flag with CCMetagen.py  __________"

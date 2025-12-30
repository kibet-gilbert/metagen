#!/bin/bash
# Script to download RefSeq FASTA files guided by ccmetagen_taxon
# Usage examples:
#   ./download_refseq.sh "bacteria"
#   ./download_refseq.sh "bacteria,fungi,viral"
#   ./download_refseq.sh "All"

set -euo pipefail

BASE_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq"
SUMMARY_FILE="assembly_summary.txt"

download_group () {
    local TAXON=$1
    local URL="${BASE_URL}/${TAXON}"
    local OUTDIR="refseq_${TAXON}_fasta"

    mkdir -p "$OUTDIR"
    cd "$OUTDIR"

    echo "Downloading assembly summary for ${TAXON}..."
    wget -q "${URL}/${SUMMARY_FILE}" -O "${SUMMARY_FILE}"

    echo "Parsing FTP links and downloading FASTA files for ${TAXON}..."
    awk -F '\t' 'BEGIN{OFS="\t"} !/^#/ {print $20}' "${SUMMARY_FILE}" | while read -r ftp_path; do
        fasta_url="${ftp_path}/${ftp_path##*/}_genomic.fna.gz"
        echo "Fetching: $fasta_url"
        wget -q "$fasta_url"
    done

    cd ..
    echo "Completed downloads for ${TAXON}"
}

# Default list of major RefSeq groups
ALL_TAXA=(archaea bacteria fungi invertebrate plant protozoa vertebrate_mammalian vertebrate_other viral)

if [ $# -eq 0 ]; then
    echo "No ccmetagen_taxon provided. Defaulting to All groups..."
    TAXA="${ALL_TAXA[*]}"
else
    TAXA="$1"
fi

# Handle "All" keyword
if [[ "$TAXA" == "All" ]]; then
    SELECTED=("${ALL_TAXA[@]}")
else
    # Split comma-separated string into array
    IFS=',' read -r -a SELECTED <<< "$TAXA"
fi

# Download each selected taxon
for TAXON in "${SELECTED[@]}"; do
    download_group "$TAXON"
done

echo "All requested RefSeq FASTA downloads complete."


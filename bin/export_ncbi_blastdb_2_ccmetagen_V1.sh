#!/usr/bin/env bash
set -euo pipefail

# # Load modules:
# module load blast/2.14.1+
# module load seqkit/2.12.0

DB_PREFIX="${1:-}"
OUTDIR="${2:-ccm_out}"

if [[ -z "$DB_PREFIX" ]]; then
  echo "Usage: $0 <blast_db_prefix> [outdir]"
  echo "Example: $0 /export/data/bio/ncbi/blast/db/v5/nt ./nt_out"
  echo "Example: $0 /export/data/ilri/sarscov/databases/ncbi/core_nt/core_nt ccm_out"
  exit 1
fi

mkdir -p "$OUTDIR"

# Direct one‑liner (skip CCMetagen’s renamer): write headers as >taxid|sequence_description straight from the BLAST DB.
echo "==> One step - Writing sequential FASTA with accession-first headers ..."
# From your BLAST db (prefix: core_nt), write desired CCMetagen-ready FASTA:
blastdbcmd -db "$DB_PREFIX" \
  -entry all \
  -outfmt ">%T|%S|%t\n%s" \
  > "${OUTDIR}/blast_nt.taxid_desc.fasta"

echo "==> Done."
echo "FASTA: ${OUTDIR}/blast_nt.taxid_desc.fasta"
# Quick sanity check: header and sequence
echo "head -n 2 ${OUTDIR}/blast_nt.taxid_desc.fasta"

# echo "==> Writing accession→taxid map ..."
# # %a = accession; %T = taxid
# # Some entries may lack an accession; keep only complete lines (2 fields)
# blastdbcmd -db "$DB_PREFIX" -entry all \
#   -outfmt "%a %T" | awk 'NF==2' > "${OUTDIR}/accession_taxid_nucl.map"
# 
# echo "==> Writing sequential FASTA with accession-first headers ..."
# # Header: >ACCESSION <TITLE>, Sequence on a single line (sequential)
# # (%a accession, %t title, %s sequence, \n newline)
# blastdbcmd -db "$DB_PREFIX" -entry all \
#   -outfmt ">%a %t\n%s" > "${OUTDIR}/nt_sequential.fa"
# 
# echo -e "==> Done.\nMap: ${OUTDIR}/accession_taxid_nucl.map"
# echo "FASTA: ${OUTDIR}/nt_sequential.fa"

# Title filter (drop “PREDICTED:”, vectors, environmental):
# Remove predicted entries:
grep -i -v '^>[^|]\+|PREDICTED:' "${OUTDIR}/blast_nt.taxid_desc.fasta" \
   | grep -i -v -E 'vector|cloning|synthetic construct|patent|uncultured|environmental sample|metagenome' \
   | awk 'BEGIN{RS=">"; ORS=""} NR>1 {print ">"$0}' \
   > "${OUTDIR}/blast_nt.taxid_desc.clean.fasta"

# Length + dedup + optional masking:
# Length filter
seqkit seq -m 300 \
   "${OUTDIR}/blast_nt.taxid_desc.clean.fasta" \
   > "${OUTDIR}/blast_nt.taxid_desc.clean.len300.fasta"
# Deduplicate
seqkit rmdup -s -i \
   "${OUTDIR}/blast_nt.taxid_desc.clean.len300.fasta" \
   > "${OUTDIR}/blast_nt.taxid_desc.clean.len300.dedub.fasta"
# optional: mask low complexity
dustmasker -in  "${OUTDIR}/blast_nt.taxid_desc.clean.len300.dedub.fasta" \
   -outfmt fasta \
   | sed 's/[a-z]/N/g' \
   > "${OUTDIR}/blast_nt.taxid_desc.clean.len300.dedub.masked.fasta"


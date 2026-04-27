#!/usr/bin/env python3
"""
Build a combined GTDB bac120 FASTA with KMA/CCMetagen-compatible headers.

Input sources (one required):
  1) NUL-delimited list of .fna file paths via STDIN
     (e.g. from: find ... -print0)
  2) --extract-dir: directory searched recursively for *.fna

Precedence:
  - If STDIN contains data, it is used
  - Otherwise --extract-dir is searched

Output FASTA header:
  >genome~marker|taxon
Usage:
========================================
  build_gtdb_scg_fasta.py \
    --extract-dir "$EXTRACT_DIR" \
    --taxonomy taxonomy_map.tsv \
    --out-fasta out.fna \
    --header-map header_map.tsv
========================== Alternatively:===========
  find "$EXTRACT_DIR" -type f -name '*.fna' -print0 \
  | build_gtdb_scg_fasta.py \
    --taxonomy taxonomy_map.tsv \
    --out-fasta out.fna \
    --header-map header_map.tsv
========================================
Author: Gilbert Kibet - GTDB SCG pipeline
"""

import sys
import os
import argparse

# ----------------------------
# Argument parsing
# ----------------------------

ap = argparse.ArgumentParser(description="Build combined GTDB bac120 FASTA")
ap.add_argument("--taxonomy", required=True, help="TSV: genome<TAB>lineage")
ap.add_argument("--out-fasta", required=True, help="Output FASTA")
ap.add_argument("--header-map", required=True, help="Header map TSV")
ap.add_argument("--extract-dir", help="Recursively search for *.fna (if no stdin)")
args = ap.parse_args()

# ----------------------------
# Load taxonomy
# ----------------------------

taxonomy = {}

if not os.path.isfile(args.taxonomy):
    sys.exit(f"ERROR: Taxonomy file not found: {args.taxonomy}")

with open(args.taxonomy) as fh:
    for line in fh:
        line = line.rstrip('\r\n')
        if not line:
            continue
        try:
            genome, lineage = line.split('\t', 1)
        except ValueError:
            sys.exit(f"ERROR: Invalid taxonomy line (expected ≥2 columns): {line}")
        taxonomy[genome] = lineage
    # Add a sanity check (do this immediately)
    print(f"Loaded {len(taxonomy)} genomes from taxonomy file", file=sys.stderr)

# ----------------------------
# Collect .fna files
# ----------------------------

fna_files = []

# Case 1: paths supplied via stdin
if not sys.stdin.isatty():
    data = sys.stdin.buffer.read()
    if data:
        fna_files = [p.decode() for p in data.split(b'\0') if p]

# Case 2: no stdin, use --extract-dir
if not fna_files and args.extract_dir:
    for root, _, files in os.walk(args.extract_dir):
        for f in files:
            if f.endswith(".fna"):
                fna_files.append(os.path.join(root, f))

# Failure case
if not fna_files:
    sys.exit(
        "ERROR: No .fna files found.\n"
        "Provide paths via stdin or use --extract-dir"
    )

fna_files.sort()

# ----------------------------
# Ensure header_map has header
# ----------------------------

HEADER_LINE = "seq_id\tgenome\tmarker\ttaxon\n"

add_header = False
if not os.path.exists(args.header_map):
    add_header = True
else:
    with open(args.header_map) as fh:
        first_line = fh.readline()
        if first_line.strip() != HEADER_LINE.strip():
            add_header = True

# ----------------------------
# Process FASTA files
# ----------------------------

try:
    fo = open(args.out_fasta, 'w')
    hm = open(args.header_map, 'a')
    if add_header:
        hm.write(HEADER_LINE)
except Exception as e:
    sys.exit(f"ERROR: Cannot open output files: {e}")

for fpath in fna_files:
    marker = os.path.splitext(os.path.basename(fpath))[0]

    with open(fpath) as fi:
        genome = None
        seq = []

        for line in fi:
            # line = line.rstrip('\r\n')
            line = line.strip()
            if line.startswith('>'):
                if genome is not None:
                    taxon = taxonomy.get(genome, "unknown")
                    seq_id = f"{genome}~{marker}"
                    fo.write(f">{seq_id}|{taxon}\n")
                    fo.write("\n".join(seq) + "\n")
                    hm.write(f"{seq_id}\t{genome}\t{marker}\t{taxon}\n")
                    if taxon == "unknown" and len(sys.argv) > 1:
                        print(f"Warning: genome '{genome}' not in taxonomy", file=sys.stderr)

                # genome = line[1:].strip()
                genome = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line)

        # Final record
        if genome is not None:
            taxon = taxonomy.get(genome, "unknown")
            seq_id = f"{genome}~{marker}"
            fo.write(f">{seq_id}|{taxon}\n")
            fo.write("\n".join(seq) + "\n")
            hm.write(f"{seq_id}\t{genome}\t{marker}\t{taxon}\n")


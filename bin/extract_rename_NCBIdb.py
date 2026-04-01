#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Add taxids to nt/core_nt/nt-like collections for CCMetagen

- Mode A: Provide --acc2tax and --infa (sequential FASTA with accession first)
- Mode B: Provide --db-prefix (BLAST DB prefix). The script will:
    1) Create accession->taxid map
    2) Create sequential FASTA with accession-first headers
    3) Then perform renaming to >taxid|sequence_description

Usage - Mode A:
python3 extract_rename_NCBIdb.py \
  --infa refseq_core_nt_out/nt_sequential.fa \
  --acc2tax refseq_core_nt_out/accession_taxid_nucl.map \
  --outfa core_nt.w_taxid.fasta

Usage - Mode B:
python3 extract_rename_NCBIdb.py \
  --db-prefix /export/data/bio/ncbi/blast/db/v5/core_nt \
  --outfa core_nt.w_taxid.fasta \
  --workdir ccm_tmp

@original: V.R. Marcelino (2018)
@updated: 2026-03-12 by G.K. with argparse, BLAST DB mode, and CCMetagen header
"""

import argparse
import logging
import os
import subprocess
import sys
from typing import Dict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

def run(cmd: list, stdout_path: str = None):
    logging.info("Running: %s", " ".join(cmd))
    if stdout_path:
        with open(stdout_path, "w") as out:
            subprocess.run(cmd, check=True, stdout=out)
    else:
        subprocess.run(cmd, check=True)

def build_from_blastdb(db_prefix: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    acc2tax = os.path.join(outdir, "accession_taxid_nucl.map")
    infa    = os.path.join(outdir, "nt_sequential.fa")

    # 1) accession -> taxid map
    # Keep only lines with two fields (accession and taxid)
    cmd_map = [
        "blastdbcmd", "-db", db_prefix, "-entry", "all",
        "-outfmt", "%a %T"
    ]
    logging.info("Creating accession→taxid map...")
    with open(acc2tax, "w") as out:
        p = subprocess.Popen(cmd_map, stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            if len(line.split()) == 2:
                out.write(line)
        p.wait()
        if p.returncode != 0:
            raise RuntimeError("blastdbcmd for map failed")

    # 2) sequential FASTA with accession-first header
    cmd_fa = [
        "blastdbcmd", "-db", db_prefix, "-entry", "all",
        "-outfmt", ">" + r"%a %t" + r"\n" + r"%s"
    ]
    logging.info("Creating sequential FASTA...")
    run(cmd_fa, stdout_path=infa)

    return acc2tax, infa

def load_map(acc2tax_path: str) -> Dict[str, str]:
    logging.info("Loading accession→taxid map from %s", acc2tax_path)
    acc2tax = {}
    with open(acc2tax_path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            acc, taxid = parts
            acc2tax[acc] = taxid
    logging.info("Loaded %d accessions", len(acc2tax))
    return acc2tax

def main():
    ap = argparse.ArgumentParser(
        description="Rename nt/core_nt FASTA headers to >taxid|sequence_description (CCMetagen style)."
    )
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--db-prefix", help="NCBI BLAST database prefix (e.g., /path/to/core_nt)")
    g.add_argument("--infa", help="Input sequential FASTA with accession-first headers")

    ap.add_argument("--acc2tax", help="accession->taxid mapping file (space-separated). Required with --infa.")
    ap.add_argument("--outfa", default="nt_w_taxid.fas", help="Output FASTA (default: nt_w_taxid.fas)")
    ap.add_argument("--workdir", default="ccm_work", help="Working directory when using --db-prefix")
    args = ap.parse_args()

    if args.db_prefix:
        # Build inputs from BLAST DB
        acc2tax_path, infa_path = build_from_blastdb(args.db_prefix, args.workdir)
    else:
        # Use provided files
        if not args.acc2tax or not args.infa:
            ap.error("When using --infa, you must also provide --acc2tax")
        acc2tax_path = args.acc2tax
        infa_path    = args.infa

    acc2tax_dic = load_map(acc2tax_path)

    logging.info("Renaming headers to >taxid|sequence_description ...")
    n_hdr = 0
    n_unmapped = 0
    with open(infa_path) as infh, open(args.outfa, "w") as outfh:
        seq_buf = []
        current_taxid = None
        current_title = None

        for line in infh:
            if line.startswith(">"):
                # flush previous
                if current_taxid is not None:
                    outfh.write(f">{current_taxid}|{current_title}\n")
                    outfh.write("".join(seq_buf))
                    seq_buf = []

                n_hdr += 1
                # Header format we created: >ACCESSION <TITLE...>
                line = line[1:].rstrip()
                if line:
                    # Split once: accession, then title (may contain spaces)
                    parts = line.split(" ", 1)
                    accession = parts[0]
                    title = parts[1] if len(parts) > 1 else ""
                    taxid = acc2tax_dic.get(accession)
                    if taxid is None:
                        taxid = "unk_taxid"
                        n_unmapped += 1
                    current_taxid = taxid
                    current_title = title
                else:
                    current_taxid = "unk_taxid"
                    current_title = ""
            else:
                # Sequence lines are already sequential (unwrapped), but keep generality
                seq_buf.append(line)

        # flush last record
        if current_taxid is not None:
            outfh.write(f">{current_taxid}|{current_title}\n")
            outfh.write("".join(seq_buf))

    logging.info("Done. Wrote %s", args.outfa)
    if n_unmapped:
        logging.warning("Unmapped accessions: %d out of %d headers", n_unmapped, n_hdr)

if __name__ == "__main__":
    main()

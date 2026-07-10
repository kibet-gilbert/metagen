#!/usr/bin/env python3
"""
extract_reads_by_taxid.py

Standalone script for per-taxid read extraction from Kraken2 results.
Stored in pipeline bin/ so Nextflow adds it to PATH automatically.

Subcommands:
  parse-report    Parse kraken2 report → taxonomy tree + abundance table
  extract-readids Filter kraken2 output → read_id/taxid pairs
  subsample       Pool read IDs per taxid group → subsample without replacement

Usage:
  extract_reads_by_taxid.py parse-report  --report <file> --taxids <csv> --prefix <str>
  extract_reads_by_taxid.py extract-readids --output <file> --taxids <file> --prefix <str>
  extract_reads_by_taxid.py subsample     --readids <file> --groups <file> --max-reads <int> --seed <int> --prefix <str>
"""

import argparse
import json
import random
import sys
from collections import defaultdict
from pathlib import Path


# =============================================================================
# SHARED UTILITIES
# =============================================================================

def info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)

def warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr)


# =============================================================================
# TAXONOMY TREE FUNCTIONS  (used by parse-report)
# =============================================================================

def parse_kraken_report(report_file):
    """
    Parse a kraken2 report file into a node dictionary.

    Kraken2 report format (tab-separated):
        pct_clade  n_reads_clade  n_reads_direct  rank_code  taxid  name

    The name column has leading spaces (2 per taxonomy depth level).
    Returns: dict mapping taxid -> node dict
    """
    nodes = {}
    stack = []   # list of (depth, taxid) to track ancestry

    with open(report_file) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue

            pct      = float(parts[0].strip())
            n_clade  = int(parts[1].strip())
            n_direct = int(parts[2].strip())
            rank     = parts[3].strip()
            taxid    = parts[4].strip()
            name_raw = parts[5]                         # has leading spaces
            depth    = (len(name_raw) - len(name_raw.lstrip(' '))) // 2
            name     = name_raw.strip()

            # Pop stack entries at same or deeper level to find true parent
            while stack and stack[-1][0] >= depth:
                stack.pop()

            parent = stack[-1][1] if stack else None

            nodes[taxid] = {
                'taxid':    taxid,
                'name':     name,
                'rank':     rank,
                'n_clade':  n_clade,
                'n_direct': n_direct,
                'pct':      pct,
                'depth':    depth,
                'parent':   parent,
                'children': [],
            }
            if parent and parent in nodes:
                nodes[parent]['children'].append(taxid)

            stack.append((depth, taxid))

    return nodes


def get_descendants(nodes, taxid):
    """
    BFS: return set of all descendant taxids for a given root taxid.
    Does not include the root taxid itself.
    """
    result = set()
    queue  = list(nodes.get(taxid, {}).get('children', []))
    while queue:
        current = queue.pop(0)
        result.add(current)
        queue.extend(nodes.get(current, {}).get('children', []))
    return result


def build_taxid_groups(nodes, targets):
    """
    For each target taxid, build a frozenset containing that taxid
    and all of its descendants.

    Returns: dict mapping target_taxid -> frozenset of all member taxids
    """
    groups = {}
    for t in targets:
        if t not in nodes:
            warn(f"taxid {t} not found in report - may have 0 reads")
            groups[t] = frozenset([t])
        else:
            desc = get_descendants(nodes, t)
            groups[t] = frozenset({t} | desc)
            info(f"taxid {t} ({nodes[t]['name']}): "
                 f"{len(desc)} descendants, "
                 f"{nodes[t]['n_clade']:,} clade reads")
    return groups


# =============================================================================
# SUBCOMMAND: parse-report
# =============================================================================

def cmd_parse_report(args):
    """
    Sub-step 1a: Parse kraken2 report.

    Outputs (all written to disk using --prefix):
      <prefix>.all_taxids.txt          flat list of all taxids of interest
      <prefix>.taxid_groups.json       target_taxid -> [member taxids]
      <prefix>.taxid_abundance.tsv     abundance table per taxid
    """
    targets = [t.strip() for t in args.taxids.split(',') if t.strip()]
    info(f"Parsing report: {args.report}")
    info(f"Target taxids: {targets}")

    nodes  = parse_kraken_report(args.report)
    groups = build_taxid_groups(nodes, targets)

    all_taxids = set()
    for s in groups.values():
        all_taxids |= s

    info(f"Total taxids of interest: {len(all_taxids)} across {len(groups)} groups")

    # 1. Flat taxid list
    all_taxids_file = f"{args.prefix}.all_taxids.txt"
    with open(all_taxids_file, 'w') as fh:
        for t in sorted(all_taxids):
            fh.write(t + '\n')

    # 2. Taxid groups JSON
    groups_file = f"{args.prefix}.taxid_groups.json"
    with open(groups_file, 'w') as fh:
        json.dump({k: sorted(v) for k, v in groups.items()}, fh, indent=2)

    # 3. Abundance table
    abundance_file = f"{args.prefix}.taxid_abundance.tsv"
    with open(abundance_file, 'w') as fh:
        cols = ['target_taxid', 'target_name', 'taxid', 'name',
                'rank', 'n_reads_clade', 'n_reads_direct', 'pct_clade']
        fh.write('\t'.join(cols) + '\n')
        for target, member_set in groups.items():
            target_name = nodes.get(target, {}).get('name', 'unknown')
            for t in sorted(member_set):
                if t in nodes:
                    n = nodes[t]
                    row = [target, target_name, t, n['name'],
                           n['rank'], str(n['n_clade']),
                           str(n['n_direct']), str(n['pct'])]
                    fh.write('\t'.join(row) + '\n')

    info(f"Written: {all_taxids_file}")
    info(f"Written: {groups_file}")
    info(f"Written: {abundance_file}")


# =============================================================================
# SUBCOMMAND: extract-readids
# =============================================================================

def cmd_extract_readids(args):
    """
    Sub-step 1b: Filter kraken2 per-read output to target taxid set.

    Kraken2 output format (tab-separated):
        C/U  read_id  taxid  length  kmer_hits

    Only classified reads (C) whose taxid is in the provided taxid list
    are written to output.

    Output: <prefix>.readid_taxid.tsv  (read_id, taxid)
    """
    # Load valid taxids
    with open(args.taxids) as fh:
        valid_taxids = set(line.strip() for line in fh if line.strip())

    info(f"Filtering kraken2 output against {len(valid_taxids):,} taxids")
    info(f"Reading: {args.output}")

    out_file = f"{args.prefix}.readid_taxid.tsv"
    n_total  = 0
    n_kept   = 0

    with open(args.output) as fin, open(out_file, 'w') as fout:
        fout.write('read_id\ttaxid\n')
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            n_total += 1
            status  = parts[0].strip()
            read_id = parts[1].strip()
            taxid   = parts[2].strip()
            if status == 'C' and taxid in valid_taxids:
                fout.write(f"{read_id}\t{taxid}\n")
                n_kept += 1

    info(f"Reads in kraken2 output  : {n_total:,}")
    info(f"Reads matching target taxids: {n_kept:,}")
    info(f"Written: {out_file}")


# =============================================================================
# SUBCOMMAND: subsample
# =============================================================================

def cmd_subsample(args):
    """
    Sub-step 1c: Pool read IDs per taxid group and subsample.

    For each target taxid group:
      - Pool all read IDs from target taxid + all descendant taxids
      - Randomly sample without replacement up to --max-reads
      - Write sampled read IDs to per_taxid/<prefix>_<taxid>.readids.txt

    The same read ID serves for both R1 and R2 (paired reads share the
    same ID in kraken2 output), so a single list is sufficient.

    Outputs:
      per_taxid/<prefix>_<target_taxid>.readids.txt  (one per group)
      <prefix>.subsampling_summary.tsv
    """
    random.seed(args.seed)
    info(f"Subsampling: max_reads={args.max_reads}, seed={args.seed}")

    # Load taxid groups
    with open(args.groups) as fh:
        taxid_groups = json.load(fh)   # {target: [member taxids]}

    # Load read_id -> taxid mapping
    taxid_to_readids = defaultdict(list)
    with open(args.readids) as fh:
        next(fh)   # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) == 2:
                read_id, taxid = parts[0], parts[1]
                taxid_to_readids[taxid].append(read_id)

    info(f"Unique taxids with reads: {len(taxid_to_readids):,}")

    # Ensure output directory exists
    Path("per_taxid").mkdir(exist_ok=True)

    summary_rows = []

    for target, member_taxids in taxid_groups.items():
        # Pool read IDs from all members of this group
        pooled = []
        for t in member_taxids:
            pooled.extend(taxid_to_readids.get(t, []))

        n_before = len(pooled)

        if n_before == 0:
            warn(f"No reads found for taxid group {target}")
            sampled = []
        elif n_before <= args.max_reads:
            sampled = pooled
        else:
            sampled = random.sample(pooled, args.max_reads)

        n_after = len(sampled)

        # Write read ID list for this group
        id_file = f"per_taxid/{args.prefix}_{target}.readids.txt"
        with open(id_file, 'w') as fh:
            for rid in sampled:
                # Write the BARE read ID — no /1 or /2 suffix.
                # The Nextflow module appends the correct suffix per
                # mate file at extraction time, since R1 and R2 need
                # different suffixes and not all FASTQ sources use them.
                fh.write(rid + '\n')

        did_subsample = 'yes' if n_before > args.max_reads else 'no'
        info(f"taxid {target}: pooled={n_before:,} → sampled={n_after:,} "
             f"(members={len(member_taxids)}, subsampled={did_subsample})")

        summary_rows.append({
            'target_taxid':     target,
            'n_member_taxids':  len(member_taxids),
            'n_reads_pooled':   n_before,
            'n_reads_sampled':  n_after,
            'subsampled':       did_subsample,
            'seed':             args.seed,
            'readid_file':      id_file,
        })

    # Write subsampling summary
    summary_file = f"{args.prefix}.subsampling_summary.tsv"
    cols = ['target_taxid', 'n_member_taxids', 'n_reads_pooled',
            'n_reads_sampled', 'subsampled', 'seed', 'readid_file']
    with open(summary_file, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for row in summary_rows:
            fh.write('\t'.join(str(row[c]) for c in cols) + '\n')

    total = sum(r['n_reads_sampled'] for r in summary_rows)
    info(f"Total reads to extract across all groups: {total:,}")
    info(f"Written: {summary_file}")


# =============================================================================
# CLI ENTRY POINT
# =============================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        prog='extract_reads_by_taxid.py',
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest='command', required=True)

    # ── parse-report ────────────────────────────────────────────────────
    p_report = sub.add_parser(
        'parse-report',
        help='Parse kraken2 report: build taxonomy tree + abundance table',
    )
    p_report.add_argument('--report',  required=True, help='kraken2 .report file')
    p_report.add_argument('--taxids',  required=True, help='Comma-separated target taxids')
    p_report.add_argument('--prefix',  required=True, help='Output file prefix')

    # ── extract-readids ─────────────────────────────────────────────────
    p_extract = sub.add_parser(
        'extract-readids',
        help='Filter kraken2 output to read_id/taxid pairs for target taxids',
    )
    p_extract.add_argument('--output', required=True, help='kraken2 .output file (per-read)')
    p_extract.add_argument('--taxids', required=True, help='File of taxids to keep (one per line)')
    p_extract.add_argument('--prefix', required=True, help='Output file prefix')

    # ── subsample ───────────────────────────────────────────────────────
    p_sub = sub.add_parser(
        'subsample',
        help='Pool read IDs per taxid group and subsample without replacement',
    )
    p_sub.add_argument('--readids',   required=True, help='readid_taxid.tsv from extract-readids')
    p_sub.add_argument('--groups',    required=True, help='taxid_groups.json from parse-report')
    p_sub.add_argument('--max-reads', required=True, type=int, help='Max reads per taxid group')
    p_sub.add_argument('--seed',      required=True, type=int, help='Random seed')
    p_sub.add_argument('--prefix',    required=True, help='Output file prefix')

    return parser


def main():
    parser = build_parser()
    args   = parser.parse_args()

    if args.command == 'parse-report':
        cmd_parse_report(args)
    elif args.command == 'extract-readids':
        cmd_extract_readids(args)
    elif args.command == 'subsample':
        cmd_subsample(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()

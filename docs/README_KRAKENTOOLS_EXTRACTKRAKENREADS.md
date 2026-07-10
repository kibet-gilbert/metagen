# `KRAKENTOOLS_EXTRACTKRAKENREADS` — Module Documentation

**Location:** `modules/local/krakentools/extractkrakenreads/main.nf`  
**Script:** `bin/extract_reads_by_taxid.py`  
**Tools:** Python 3 · seqtk · pigz

---

## Purpose

Extracts and subsamples reads from Kraken2 classification results for a list of target pathogen taxids. Unlike `extract_kraken_reads.py` (which applies a single global read cap across all taxids combined), this module processes each target taxid group independently, ensuring that an abundant organism cannot exhaust the read budget before rarer target taxa are sampled.

---

## Sub-steps

### Step 1a — Parse kraken2 report

**Script subcommand:** `extract_reads_by_taxid.py parse-report`

Parses `*.kraken2.report.txt` to build an in-memory taxonomy tree using the indentation depth of the name column (2 spaces per level). For each target taxid, all descendant taxids are collected via breadth-first search. Outputs:

| File | Description |
|------|-------------|
| `<prefix>.all_taxids.txt` | Flat list of all taxids of interest (targets + all descendants) |
| `<prefix>.taxid_groups.json` | Mapping of each target taxid to its full member set |
| `<prefix>.taxid_abundance.tsv` | Per-taxid read counts and percentages from the report |

---

### Step 1b — Extract read IDs from kraken2 output

**Script subcommand:** `extract_reads_by_taxid.py extract-readids`

Scans `*.kraken2.output.txt` (the per-read classification file) and retains only classified reads (`C` in column 1) whose assigned taxid appears in `all_taxids.txt`. Output:

| File | Description |
|------|-------------|
| `<prefix>.readid_taxid.tsv` | Two-column table: `read_id` and assigned `taxid` |

This step operates on read IDs only — no FASTQ access yet — making it fast regardless of file size.

---

### Step 1c — Subsample read IDs per taxid group

**Script subcommand:** `extract_reads_by_taxid.py subsample`

For each target taxid group, all read IDs from that taxid and its descendants are pooled. If the pool exceeds `max_reads_per_taxid`, `random.sample()` draws without replacement using a fixed seed for reproducibility. Because Kraken2 assigns the same ID to both paired reads, a single ID list covers R1 and R2 — mate pairing is preserved by design. Outputs:

| File | Description |
|------|-------------|
| `per_taxid/<prefix>_<taxid>.readids.txt` | Sampled read ID list, one file per target group |
| `<prefix>.subsampling_summary.tsv` | Before/after read counts and subsample status per group |

---

### Step 1d — Physical read extraction

**Tool:** `seqtk subseq`

For each per-taxid read ID list, `seqtk subseq` extracts matching reads from the R1 (and R2) FASTQ files. Outputs are appended into combined FASTA files and compressed with `pigz`.

| File | Description |
|------|-------------|
| `<prefix>.extracted_R1.fasta.gz` | All extracted R1 reads (all taxid groups combined) |
| `<prefix>.extracted_R2.fasta.gz` | Matching R2 reads (paired, same IDs as R1) |

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_reads_per_taxid` | `50000` | Maximum reads extracted per target taxid group |
| `subsample_seed` | `42` | Random seed for reproducible subsampling |
| `include_children` | `true` | Include descendant taxids in each group |
| `fastq_output` | `false` | Output FASTQ instead of FASTA |

Set via `ext.*` in `nextflow.config`:

```groovy
withName: 'KRAKENTOOLS_EXTRACTKRAKENREADS' {
    ext.max_reads_per_taxid = 50000
    ext.subsample_seed      = 42
    ext.include_children    = true
    ext.fastq_output        = false
    container = "/path/to/krakentools_apptainer.sif"
}
```

---

## Output summary

```
<prefix>.extracted_R1.fasta.gz       ← BLAST input (Step 3)
<prefix>.extracted_R2.fasta.gz       ← paired R2 reads
<prefix>.taxid_abundance.tsv         ← read counts per taxid from report
<prefix>.subsampling_summary.tsv     ← before/after counts per group
<prefix>.readid_taxid.tsv            ← read_id → taxid mapping (audit trail)
<prefix>.taxid_groups.json           ← target → descendant taxid mapping
<prefix>.summary.txt                 → consumed by VALIDATE_HITS
```

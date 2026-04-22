# Building CCMetagen / KMA Databases from NCBI Nucleotide Sources

---

## Table of Contents

- [Overview](#overview)
- [Database sources](#database-sources)
- [Prerequisites](#prerequisites)
- [Part A — NCBI BLAST databases (nt and core\_nt)](#part-a--ncbi-blast-databases-nt-and-core_nt)
  - [A1. Download BLAST database](#a1-download-blast-database)
  - [A2. Export raw FASTA — Module 1: export\_blastdb.sh](#a2-export-raw-fasta--module-1-export_blastdbsh)
  - [A3. Clean FASTA — Module 2: clean\_fasta.sh](#a3-clean-fasta--module-2-clean_fastash)
- [Part B — NCBI RefSeq reference genomes](#part-b--ncbi-refseq-reference-genomes)
  - [B1–B12. Full pipeline: download\_all\_refseq.sh](#b1b12-full-pipeline-download_all_refseqsh)
  - [B13. Optional cleaning pass with clean\_fasta.sh](#b13-optional-cleaning-pass-with-clean_fastash)
- [Part C — KMA indexing](#part-c--kma-indexing)
- [Part D — NCBI database selection guide](#part-d--ncbi-database-selection-guide)
- [Troubleshooting](#troubleshooting)
- [File inventory](#file-inventory)

---

## Overview

[CCMetagen](https://github.com/vrmarcelino/CCMetagen) uses [KMA](https://bitbucket.org/genomicepidemiology/kma) for competitive read mapping against a reference database. KMA requires the reference sequences to be:

1. In FASTA format with **CCMetagen-compatible headers**: `>taxid|title`
2. Indexed with `kma index`

Three database sources are supported in this workflow:

| Source | Script(s) | When to use |
|--------|-----------|-------------|
| NCBI BLAST `core_nt` | `export_blastdb.sh` → `clean_fasta.sh` | Broad metagenomic profiling; balanced size vs sensitivity |
| NCBI BLAST `nt` | `export_blastdb.sh` → `clean_fasta.sh` | Maximum sensitivity; very large (>300 GB) |
| NCBI RefSeq reference genomes | `download_all_refseq.sh` → *(optional)* `clean_fasta.sh` | High-quality, curated genomes only; best for known pathogens |

The two-module design (`export_blastdb.sh` + `clean_fasta.sh`) means the cleaning step can be applied to **any** FASTA source — both BLAST exports and the merged RefSeq output.

---

## Database sources

### NCBI BLAST: `nt` vs `core_nt`

| Property | `nt` | `core_nt` |
|----------|------|-----------|
| Size (compressed) | ~1300 GB | ~350 GB |
| Sequences | ~100 million | ~30 million |
| Coverage | All GenBank, RefSeq, EMBL, DDBJ, patents, environmental | Curated high-quality subset |
| Redundancy | High | Reduced |
| Best for | Maximum sensitivity | Speed + quality balance |

### NCBI RefSeq (via `download_all_refseq.sh`)

RefSeq reference genomes are manually curated, one per species or strain. They are the highest-quality sequences in the NCBI ecosystem but cover fewer organisms than `nt` or `core_nt`. For wastewater metagenomics, RefSeq is ideal for pathogen detection while `core_nt` or `nt` captures the broader microbial community.

---

## Prerequisites

### Software

```bash
# NCBI BLAST+ (for blastdbcmd and update_blastdb.pl)
module load blast/2.14.1+   # or equivalent on your HPC

# KMA (for indexing)
# https://bitbucket.org/genomicepidemiology/kma
module load kma

# seqkit (for length filtering and deduplication)
# https://bioinf.shenwei.me/seqkit/
conda install -c bioconda seqkit

# GNU Parallel
conda install -c conda-forge parallel

# pigz (parallel gzip; falls back to gzip if absent)
conda install -c conda-forge pigz

# dustmasker (optional; for low-complexity masking; part of BLAST+)
# included with blast/2.14.1+

# python3 (for eutils taxid lookup in download_all_refseq.sh)
module load python/3.9
```

### Disk space estimates

| Database | Raw export | After cleaning | KMA index |
|----------|------------|---------------|-----------|
| `core_nt` | ~400 GB | ~360–380 GB | +200–300 GB |
| `nt` | ~1300 GB | ~1180–1220 GB | +600 GB |
| RefSeq reference | ~350 GB (raw .fna.gz) | ~200 GB (merged) | +50 GB |

> **Tip:** Run on a scratch filesystem with at least 2× the final size free, to accommodate intermediates.

---

## Part A — NCBI BLAST databases (nt and core\_nt)

### A1. Download BLAST database

BLAST databases are downloaded using `update_blastdb.pl`, which is bundled with BLAST+.

```bash
# List all available databases
module load blast/2.14.1+
update_blastdb.pl --showall

# Download core_nt (recommended starting point)
mkdir -p /export/data/ncbi/blast/db/v5/core_nt
cd       /export/data/ncbi/blast/db/v5/core_nt
update_blastdb.pl --decompress core_nt
# Runtime: 2–6 hours depending on network speed

# Or download full nt (larger)
mkdir -p /export/data/ncbi/blast/db/v5/nt
cd       /export/data/ncbi/blast/db/v5/nt
update_blastdb.pl --decompress nt
# Runtime: 6–24 hours

# Verify download
blastdbcmd -db core_nt -info
```

The result is a directory of `.nsq`, `.nin`, `.ndb`, `.nhr` files (BLAST v5 format) plus taxonomy metadata (`taxdb.bti`, `taxdb.btd`).

---

### A2. Export raw FASTA — Module 1: `export_blastdb.sh`

This script exports all sequences from the BLAST database with CCMetagen-compatible headers (`>taxid|title`) and optionally compresses the output.

```bash
# Syntax
export_blastdb.sh <blast_db_prefix> [outdir]

# Export core_nt
export_blastdb.sh \
  /export/data/ncbi/blast/db/v5/core_nt/core_nt \
  ./core_nt_export

# Export full nt
export_blastdb.sh \
  /export/data/ncbi/blast/db/v5/nt/nt \
  ./nt_export

# Override header format (include scientific name):
HEADER_FMT=">%T|%S|%t\n%s" \
  export_blastdb.sh /path/to/nt ./nt_export
```

**Header format tokens:**

| Token | Meaning |
|-------|---------|
| `%T` | Taxid |
| `%t` | Sequence title (description line) |
| `%S` | Scientific name |
| `%s` | Sequence (nucleotides) |

Default format `>%T|%t\n%s` produces:
```
>9606|Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly
ATCGATCG...
```

**Output:**
```
core_nt_export/
├── core_nt.raw.fasta.gz    ← pass this to clean_fasta.sh
└── export_blastdb.log
```

**Runtime:** 2–8 hours for `core_nt`; 6–20 hours for `nt`.

---

### A3. Clean FASTA — Module 2: `clean_fasta.sh`

This script takes the raw exported FASTA and applies four sequential cleaning steps. It accepts both plain and gzip-compressed input, and works on BLAST exports **or** the merged RefSeq output from `download_all_refseq.sh`.

```bash
# Syntax
clean_fasta.sh <input.fasta[.gz]> [outdir]

# Clean core_nt export
clean_fasta.sh \
  ./core_nt_export/core_nt.raw.fasta.gz \
  ./core_nt_clean

# Clean nt export with custom settings
MIN_LEN=200 MASK_LOW_COMPLEX=0 THREADS=8 \
  clean_fasta.sh \
  ./nt_export/nt.raw.fasta.gz \
  ./nt_clean

# Clean RefSeq merged FASTA (see Part B)
clean_fasta.sh \
  ./refseq/refseq_ccmetagen.fna.gz \
  ./refseq_clean
```

#### What each cleaning step removes

**Step 1 — Title filtering**

Removes sequences whose header line matches any of:

| Pattern | Reason to exclude |
|---------|-------------------|
| `PREDICTED:` | Computational predictions; redundant isoforms, not experimental |
| `vector` | Cloning vectors; cause false-positive hits |
| `synthetic construct` | Lab constructs; not from biological organisms |
| `cloning` | Cloning artifacts |
| `patent` | Patent sequences; often redundant synthetic sequences |
| `uncultured` | Unclassified environmental sequences; noisy taxonomy |
| `environmental sample` | Poor taxonomic resolution |
| `metagenome` | Circular reference — metagenomic sequences in a metagenomics DB |

**Step 2 — Length filter (default ≥ 300 nt)**

Short sequences (<300 nt) produce unreliable taxonomic assignments and inflate multi-mapper counts. Adjust with `MIN_LEN=200` if you need to capture shorter amplicons or viral genomes.

**Step 3 — Deduplication (`seqkit rmdup -s -i`)**

Removes sequences with identical nucleotide content regardless of header. Identical sequences from different submitters (common in nt) inflate database size without adding information. A list of removed duplicates is written to `duplicates.txt` for audit.

**Step 4 — Low-complexity masking (`dustmasker`, optional)**

Masks repetitive / low-complexity regions (poly-A tails, simple repeats, satellite sequences) to uppercase `N`. This prevents these regions from driving alignment scores and causing spurious hits. Disable with `MASK_LOW_COMPLEX=0` if you need exact sequences preserved.

**Output:**
```
core_nt_clean/
├── core_nt.raw.ccmetagen.masked.fasta.gz   ← ready for kma index
├── duplicates.txt
├── duplicate_seqids.txt
└── clean_fasta.log
```

**Environment variables for `clean_fasta.sh`:**

| Variable | Default | Description |
|----------|---------|-------------|
| `MIN_LEN` | 300 | Minimum sequence length (nt) |
| `MASK_LOW_COMPLEX` | 1 | 1=run dustmasker, 0=skip |
| `COMPRESS_FINAL` | 1 | 1=gzip output, 0=plain FASTA |
| `KEEP_INTERMEDIATE` | 0 | 1=keep step-by-step tmp files |
| `THREADS` | 4 | Threads for seqkit/pigz |

---

## Part B — NCBI RefSeq reference genomes

RefSeq reference genomes are downloaded directly from NCBI FTP rather than from the BLAST database. The header format is stamped post-download to match CCMetagen requirements.

### B1–B12. Full pipeline: `download_all_refseq.sh`

The script runs 12 sequential steps. All steps are idempotent — re-running the script will skip already-completed steps based on file existence checks.

```bash
# Syntax
download_all_refseq.sh [workdir] [nthreads] [want_rep]

# Download reference genomes only (default)
download_all_refseq.sh refseq 20 0

# Include representative genomes as well
download_all_refseq.sh refseq 20 1

# Keep intermediate files (useful for debugging)
KEEP_INTERMEDIATES=1 download_all_refseq.sh refseq 20 0

# Skip the final merge (if disk space is tight)
SKIP_MERGE=1 download_all_refseq.sh refseq 20 0
```

**Environment variables:**

| Variable | Default | Description |
|----------|---------|-------------|
| `KEEP_INTERMEDIATES` | 0 | 1=keep .fna.gz and stamped files after merge |
| `SKIP_MERGE` | 0 | 1=do not create merged FASTA |
| `NCBI_EMAIL` | `$USER@ilri.org` | Email for polite NCBI eutils requests |

#### Step-by-step breakdown

**Step 1 — Download assembly summary**

Downloads `assembly_summary_refseq.txt` from NCBI FTP. This file lists every assembly in RefSeq with its category, FTP path, and metadata. The full file has ~2 million entries.

```
Output: refseq/assembly_summary.txt
```

**Step 2 — Filter to reference genome URLs**

Filters the assembly summary to retain only `reference genome` entries (the highest quality, manually curated assemblies). With `want_rep=1`, `representative genome` entries are also included.

Column 5 (`refseq_category`) is filtered; column 20 (`ftp_path`) is extracted.

```
Output: refseq/assembly_summary_ref.txt
# Typical counts:
#   reference genome only:   ~25,000 URLs
#   + representative genome: ~300,000 URLs
```

**Step 3 — Parallel genome download with retries**

Downloads each genome as `{basename}_genomic.fna.gz` using `GNU Parallel` with:
- Up to `$nthreads` concurrent jobs
- Polite 0.25s delay between job starts
- Up to 10 retry rounds with exponential backoff (5s, 10s, 20s... capped at 300s)
- Per-round MD5 checksum recording
- Skip logic for already-downloaded files

Each retry round processes only the failed URLs from the previous round, so transient network errors (HTTP 429, 503) are automatically recovered.

```
Output:
  refseq/{GCF_*}_genomic.fna.gz          (one per genome)
  refseq/_dl/urls.r*.txt                 (retry lists)
  refseq/_dl/wget.r*.log                 (per-round wget log)
  refseq/_dl/parallel.r*.joblog          (GNU Parallel job log)
  refseq/_dl/md5checksum.r*.txt          (MD5 checksums)
  refseq/download_failures.final.txt     (if persistent failures remain)
```

**Step 4 — Download accession-to-taxid maps**

Downloads two NCBI taxonomy mapping files:

```
nucl_gb.accession2taxid.gz   — GenBank accessions
nucl_wgs.accession2taxid.gz  — WGS project accessions
```

These files map every NCBI nucleotide accession to its NCBI taxid. Combined they cover ~137 million accession-taxid pairs.

```bash
# Build combined versioned_accession → taxid map
zgrep "_" refseq/nucl_*.accession2taxid.gz \
  | cut -f2-3 > refseq/accession_taxid_nucl-rs.map
# ~137 million entries
```

**Step 5 — Build accession → filepath map**

Reads the header (`>` lines) of every downloaded genome FASTA to extract accession IDs, then maps each accession to its source file path. Multiple accessions per file are common (multi-chromosome genomes, plasmids).

```
Output: refseq/accession_fqpath-rs.map
# Example entry:
NC_000001.11    refseq/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
NC_000002.12    refseq/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

**Step 6 — Join maps → taxid\_fqpath.map**

Performs a sort-merge join between the accession-taxid map and the accession-filepath map on the accession column. This produces the critical mapping file: `taxid → filepath`.

```
Output: refseq/taxid_fqpath.map
# Example entries:
9606    refseq/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
287     refseq/GCF_000006765.1_ASM676v1_genomic.fna.gz
# ~25,000–26,000 entries expected
```

> **Why the join can produce fewer entries than downloaded files:**
> Some accessions in the genome FASTA headers are not present in `nucl_gb` or `nucl_wgs` accession2taxid files. This affects primarily very new RefSeq entries (submitted after the last accession2taxid release) and some non-standard accession prefixes. These are resolved in Step 9.

**Step 7 — Stamp taxID into FASTA headers**

Rewrites every FASTA header from:
```
>NC_000001.11 Homo sapiens chromosome 1, GRCh38
```
to CCMetagen-compatible format:
```
>9606|NC_000001.11 Homo sapiens chromosome 1, GRCh38
```

Uses `gzip -dc | sed | pigz -p 2 -c` in a single streaming pass — no intermediate uncompressed files. Runs in parallel with one `pigz` instance per genome using 2 threads each.

Output files: `{GCF_*}_genomic_taxID2.fna.gz`

**Step 8 — Audit: cross-reference files vs maps**

Compares three sets:
1. `.fna.gz` files on disk
2. Entries in `taxid_fqpath.map`
3. Stamped `_taxID2.fna.gz` outputs

Reports two anomaly classes:

| Code | Meaning | Action |
|------|---------|--------|
| `MAP_NO_FILE` | In map but not on disk | File download failed — re-run Step 3 |
| `FILE_NO_MAP` | On disk but not in map | Accession missing from taxid mapping — resolved in Step 9 |

```
Output: refseq/audit.txt
```

**Step 9 — Resolve unmapped genomes via NCBI eutils**

For each `FILE_NO_MAP` genome (files on disk with no taxid in `taxid_fqpath.map`), the script:

1. Extracts the GCF accession from the filename
2. Queries `esearch.fcgi` to get the NCBI Assembly database UID
3. Queries `esummary.fcgi` to retrieve the taxid
4. Appends resolved entries to `taxid_fqpath.residual_fs.map`
5. Stamps those genomes using the same single-pass method as Step 7

A 0.4-second delay between eutils requests respects NCBI's rate limit (max 3 requests/second without API key; add `&api_key=YOUR_KEY` for 10/second).

```
Output:
  refseq/missing_gcf_accessions.txt   (GCF accessions that were unmapped)
  refseq/missing_taxids.txt           (resolved taxid per GCF, or NA)
  refseq/taxid_fqpath.residual_fs.map (supplementary taxid→path map)
```

**Step 10 — Final stamp audit**

Counts input `.fna.gz` files, stamped `_taxID2.fna.gz` files, and identifies any remaining unstamped genomes. The unstamped list is written to `refseq/unstamped.txt` for manual investigation.

```
Expected output (reference genomes only):
  Input .fna.gz:            25,265
  Stamped _taxID2.fna.gz:   25,265   ← ideally 100%
  Missing:                       0
```

**Step 11 — Merge all stamped FASTA into one file**

Concatenates all `*_taxID2.fna.gz` files into a single compressed FASTA. This is the file that will be indexed by KMA.

```bash
cat refseq/*_taxID2.fna.gz > refseq/refseq_ccmetagen.fna.gz
```

```
Output: refseq/refseq_ccmetagen.fna.gz
# Expected: ~25,000 genomes, >11 million sequences
# Size: ~200–350 GB
```

**Step 12 — Cleanup intermediates**

If `KEEP_INTERMEDIATES=0` (the default), removes:
- Individual `*_genomic.fna.gz` source downloads
- Individual `*_taxID2.fna.gz` stamped files (already merged)
- Per-round download logs in `refseq/_dl/`

Retained permanently:
```
refseq/refseq_ccmetagen.fna.gz          ← merged FASTA for KMA
refseq/taxid_fqpath.map                 ← primary taxid map
refseq/taxid_fqpath.residual_fs.map     ← supplementary taxid map
refseq/accession_taxid_nucl-rs.map      ← raw accession-taxid map
refseq/assembly_summary.txt             ← assembly metadata
refseq/assembly_summary_ref.txt         ← filtered URL list
refseq/audit.txt                        ← cross-reference audit
refseq/download_failures.final.txt      ← failed downloads (if any)
refseq/missing_*.txt                    ← unmapped genome records
refseq/download_refseq.log              ← full run log
```

> **Safety guard:** The script checks that `refseq_ccmetagen.fna.gz` exists and is non-empty before removing intermediate files. If the merge failed for any reason, intermediates are preserved.

---

### B13. Optional cleaning pass with `clean_fasta.sh`

The merged RefSeq FASTA from Step 11 already has correct CCMetagen headers and consists of curated reference genomes, so aggressive filtering is less critical than for `nt`. However, `clean_fasta.sh` can still be applied to:

- Remove any remaining predicted/environmental entries in edge cases
- Apply length filtering (removes very short contigs/fragments)
- Deduplicate sequences that are identical across different assemblies

```bash
# Apply clean_fasta.sh to merged RefSeq output
MIN_LEN=200 MASK_LOW_COMPLEX=0 \
  clean_fasta.sh \
  ./refseq/refseq_ccmetagen.fna.gz \
  ./refseq_clean
```

For RefSeq, `MASK_LOW_COMPLEX=0` is recommended since curated reference genomes have already been processed by NCBI and low-complexity masking can interfere with small viral or bacterial genomes.

---

## Part C — KMA indexing

Once any of the three cleaned FASTAs are ready, index them for KMA:

```bash
# Load KMA
module load kma   # or: conda activate ccmetagen_env

# Index core_nt
kma index \
  -i ./core_nt_clean/core_nt.raw.ccmetagen.masked.fasta.gz \
  -o /export/data/kma_db/core_nt_ccmetagen \
  -verbose

# Index RefSeq
kma index \
  -i ./refseq/refseq_ccmetagen.fna.gz \
  -o /export/data/kma_db/refseq_ccmetagen \
  -verbose

# Index nt (large — may take 12–24h and require >100 GB RAM)
kma index \
  -i ./nt_clean/nt.raw.ccmetagen.masked.fasta.gz \
  -o /export/data/kma_db/nt_ccmetagen \
  -verbose
```

The index consists of several files with the given prefix:
```
refseq_ccmetagen.index    refseq_ccmetagen.seq
refseq_ccmetagen.name     refseq_ccmetagen.comp.b
refseq_ccmetagen.length.b refseq_ccmetagen.decon.comp.b
```

**Run CCMetagen with the indexed database:**

```bash
# KMA alignment
kma \
  -ipe sample_R1.fastq.gz sample_R2.fastq.gz \
  -o sample_kma_output \
  -t_db /export/data/kma_db/refseq_ccmetagen \
  -1t1 -mem_mode -and \
  -t 16

# CCMetagen taxonomic classification
CCMetagen.py \
  -i sample_kma_output.res \
  -r nt \              # use -r nt for RefSeq and nt/core_nt databases
  -o sample_ccmetagen \
  --depth 0.2
```

---

## Part D — NCBI database selection guide

| Scenario | Recommended database | Reason |
|----------|---------------------|--------|
| Broad wastewater metagenomics — all kingdoms | `core_nt` | Good balance of sensitivity and speed; includes bacteria, viruses, fungi, eukaryotes |
| Pathogen-focused surveillance | RefSeq reference genomes | Highest quality; one genome per pathogen; minimal noise |
| Maximum sensitivity (novel organisms) | `nt` | Covers all submitted sequences; very slow |
| Bacteria + archaea only | `nt_prok` | Smaller; focused on prokaryotes |
| Viral surveillance | `nt_viruses` or `ref_viruses_rep_genomes` | Viral-only databases |
| Rapid classification, less RAM | RefSeq reference | Smaller index; faster |

### Combining databases

For wastewater AMR surveillance, the recommended strategy is to build **two separate KMA databases** and run them independently:

```
1. RefSeq reference genomes  → pathogen detection (high confidence)
2. core_nt                   → community profiling (broad coverage)
```

This avoids the need to build and maintain the full `nt` database while covering both use cases.

---

## Troubleshooting

### Download failures in `download_all_refseq.sh`

```bash
# Check which URLs failed
cat refseq/download_failures.final.txt

# Re-run the script — it will only retry failed URLs
download_all_refseq.sh refseq 20 0

# If persistent 429 errors: reduce parallel jobs
download_all_refseq.sh refseq 4 0
```

### Audit shows FILE_NO_MAP entries

These are genomes downloaded successfully but with accessions not found in `nucl_gb` or `nucl_wgs` maps. Step 9 resolves most of these via eutils. If some remain as `NA` in `missing_taxids.txt`:

```bash
# Check what accession prefixes they have
grep NA refseq/missing_taxids.txt | while read _ gcf; do
  ls refseq/${gcf}_* 2>/dev/null
done | head -5 | xargs -I{} zcat {} | grep '^>' | head -3
```

New `NZ_` or `GCF_` accessions may require the `dead_nucl.accession2taxid.gz` map:
```bash
wget -O refseq/dead_nucl.accession2taxid.gz \
  https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz
```

### `blastdbcmd` export is very slow

`blastdbcmd -entry all` is sequential by design. For `nt`, expect 8–24 hours. Parallelisation is not straightforward because the BLAST DB format does not support random-access by index in a portable way. Use `core_nt` to save time.

### KMA indexing runs out of memory

```bash
# Check available memory
free -h

# For very large databases, use batch indexing
kma index -i input.fasta.gz -o db_prefix -batch 10000000
```

### `seqkit rmdup` is slow on very large files

For `nt`-scale databases (>100 million sequences), exact sequence deduplication with `seqkit rmdup -s` requires loading all sequences into memory. As an alternative, deduplicate by sequence ID only:

```bash
seqkit rmdup -j 8 -n input.fasta > dedup.fasta   # by ID only, faster
```

Or use `vsearch` for clustering at 95% identity (reduces size more aggressively):

```bash
vsearch --cluster_fast input.fasta \
  --id 0.95 \
  --centroids dedup_95.fasta \
  --threads 16
```

---

## File inventory

### Scripts

| Script | Module | Input | Output |
|--------|--------|-------|--------|
| `export_blastdb.sh` | 1 of 2 (BLAST path) | BLAST DB prefix | `*.raw.fasta.gz` |
| `clean_fasta.sh` | 2 of 2 (both paths) | any `.fasta[.gz]` | `*.ccmetagen[.masked].fasta.gz` |
| `download_all_refseq.sh` | 1+2 (RefSeq path) | — | `refseq_ccmetagen.fna.gz` |

### Key output files

| File | Description |
|------|-------------|
| `{db}.raw.fasta.gz` | Raw export from BLAST DB (before cleaning) |
| `{db}.ccmetagen.masked.fasta.gz` | Cleaned, masked, compressed FASTA ready for KMA |
| `refseq/refseq_ccmetagen.fna.gz` | Merged, header-stamped RefSeq FASTA |
| `refseq/taxid_fqpath.map` | Primary taxid → filepath map |
| `refseq/taxid_fqpath.residual_fs.map` | Supplementary taxid map (eutils-resolved) |
| `refseq/audit.txt` | Cross-reference audit of files vs maps |
| `refseq/download_failures.final.txt` | Persistent download failures (if any) |
| `refseq/missing_taxids.txt` | GCF → taxid lookups from eutils |
| `refseq/download_refseq.log` | Full run log with timestamps |
| `{outdir}/clean_fasta.log` | Step-by-step cleaning log |
| `{outdir}/duplicates.txt` | List of removed duplicate sequences |
| `{kma_db_prefix}.*` | KMA index files ready for CCMetagen |

---

*For CCMetagen usage after database building, refer to the [CCMetagen documentation](https://github.com/vrmarcelino/CCMetagen).*  
*For KMA aligner options, refer to the [KMA documentation](https://bitbucket.org/genomicepidemiology/kma).*

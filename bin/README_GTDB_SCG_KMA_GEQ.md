# GTDB bac120 Single-Copy Genes: KMA Database & Genome Equivalent Normalisation

---

## Table of Contents

- [1. Understanding your GTDB downloads](#1-understanding-your-gtdb-downloads)
- [2. Which sequence file to use](#2-which-sequence-file-to-use)
- [3. GEQ normalisation theory](#3-geq-normalisation-theory)
- [4. Building the KMA database — step by step](#4-building-the-kma-database--step-by-step)
- [5. KMA alignment strategy](#5-kma-alignment-strategy)
- [6. Interpreting KMA output for GEQ](#6-interpreting-kma-output-for-geq)
- [7. Per-SCG abundance and taxonomic abundance](#7-per-scg-abundance-and-taxonomic-abundance)
- [8. Aligning with ResFinder output format](#8-aligning-with-resfinder-output-format)
- [9. References](#9-references)

---

## 1. Understanding your GTDB downloads

### `bac120_msa_marker_info_r232.tsv` — Marker gene metadata

This is a plain TSV describing all 120 bacterial single-copy marker genes (SCGs) used to build the GTDB phylogenetic tree. Each row is one marker.

```
Marker Id         Name          Description                           Length  Single copy (%)  Ubiquity (%)
PFAM_PF00380.20   Ribosomal_S9  Ribosomal protein S9/S16             121     91.47            92.77
TIGR_TIGR00006    TIGR00006     16S rRNA methyltransferase            310     90.94            95.08
```

The columns explain:

| Column | Meaning | Relevance for GEQ |
|--------|---------|-------------------|
| `Marker Id` | Pfam or TIGRFAM accession + version | Unique database ID |
| `Name` | Short gene name | Used to label output |
| `Description` | Full gene function | Annotation |
| `Length (bp)` | Median alignment length in the MSA | Used to normalise RPKM denominator |
| `Single copy (%)` | % of GTDB genomes where this gene appears exactly once | **Key quality filter** — only use genes ≥85% |
| `Ubiquity (%)` | % of GTDB genomes where this gene appears at all | **Key quality filter** — only use genes ≥85% |

For GEQ normalisation, you want markers that are:
- Present in **nearly all bacteria** (high ubiquity) — otherwise the estimate is biased towards taxa that have those markers
- Present in **exactly one copy per genome** (high single copy %) — so that read depth reflects genome count, not gene copy number

**Recommended filter:** `Single copy (%) >= 85 AND Ubiquity (%) >= 85`. From the preview above this retains most of the 120 markers. A stricter filter of ≥90% for both retains approximately 40–60 markers and gives the most reliable GEQ estimates.

---

### `bac120_metadata_r232.tsv.gz` — Genome-level metadata

A very wide TSV (~100 columns) with one row per genome in GTDB r232. Key columns:

| Column | Meaning |
|--------|---------|
| `accession` | GTDB accession (e.g. `RS_GCF_027889305.1`) |
| `gtdb_taxonomy` | Full 7-rank GTDB taxonomy string |
| `gtdb_genome_representative` | Accession of the species cluster representative |
| `gtdb_representative` | `t/f` — is this genome the cluster representative? |
| `checkm2_completeness` | Genome completeness (%) from CheckM2 |
| `checkm2_contamination` | Estimated contamination (%) |
| `genome_size` | Total assembly size (bp) |
| `ncbi_organism_name` | NCBI organism name |
| `ncbi_refseq_category` | reference/representative/na |

This file is used to:
1. Map KMA hits back to GTDB taxonomy
2. Filter genomes by quality (`checkm2_completeness ≥ 95, contamination ≤ 5`)
3. Identify which genomes are representative (one per species cluster)

---

### `bac120_marker_genes_reps_r232.tar.gz` — **Recommended for GEQ**

Contains nucleotide sequences of all 120 bac120 marker genes extracted from **representative genomes only** — one genome per GTDB species cluster (~80,000 species in r232).

Archive structure after extraction:
```
bac120_marker_genes_reps_r232/
├── faa                                                                         
│   ├── PF00380.20.faa      ← one FASTA per marker (multiple genomes)
│   ├── PF00410.20.faa
│   ├── ...
│   ├── TIGR03725.faa
│   └── ... (120 files, one per marker)
└── fna
    ├── PF00380.20.fna
    ├── PF00410.20.fna
    ├── ...
    ├── TIGR03725.fna
    └── ... (120 files, one per marker)
```

Each file is given the `marker ID` name: `*.faa` - for amino acid sequences and `*.fna` for nuclueotide sequences.  
Each FASTA header follows the format:
```
>RS_GCF_000001405.40
[gene_start]-[gene_end]
```

- `RS_GCF_000001405.40` — GTDB genome accession
- Optional: `*_1` — copy number index (copies >1 indicate multi-copy instances to filter)
- Filename: `PF00380.20` — marker ID

**Size:** ~9–15 GB compressed; ~30–50 GB extracted.

---

### `bac120_marker_genes_all_r232.tar.gz` — All genomes (not recommended for GEQ)

Same structure as above but includes **all ~700,000 genomes** in GTDB r232, not just representatives. This is much larger (~15 GB compressed) and contains near-identical sequences for closely related strains.

**Do not use this for GEQ** — the redundant sequences will cause KMA to assign reads to multiple near-identical templates, inflating multi-mapper counts and distorting abundance estimates. Use `_reps_` only.

---

## 2. Which sequence file to use

**Use `bac120_marker_genes_reps_r232.tar.gz`** for KMA-based GEQ normalisation.

| File | Genomes | Use case |
|------|---------|----------|
| `_reps_` | ~80,000 representative genomes | **GEQ normalisation** ✓ |
| `_all_` | ~700,000 all genomes | Strain-level diversity studies only |

The representative set ensures each species contributes one sequence per marker. KMA's competitive mapping will assign reads to the best-matching species, and the abundance per marker directly reflects genome copies.

---

## 3. GEQ normalisation theory

### Why RPKM is insufficient for cross-sample ARG comparisons

RPKM normalises for sequencing depth (total reads) and gene length, but does not correct for variation in **microbial biomass** between samples. Two samples with the same RPKM for a resistance gene could have very different biological meanings:

- Sample A: 1,000 RPKM ARG, 10 genome equivalents per ml → 100 ARG copies per genome
- Sample B: 1,000 RPKM ARG, 100 genome equivalents per ml → 10 ARG copies per genome

### Genome Equivalent (GEQ) definition

Universal single-copy genes (USiCGs) are found in the genomes of almost all microorganisms, generally at a single copy, and therefore represent a proxy for invariant genomic elements. The GTDB bac120 markers are precisely such genes.

A **genome equivalent** is defined as:

```
GEQ = median(RPKM_per_SCG_marker)  across all qualifying SCG markers in a sample
```

This estimates how many bacterial genome copies are represented per million reads, normalising for:
- Sequencing depth (via RPKM)
- Gene length (via kilobase normalisation)
- The fact that one genome contributes exactly one copy of each SCG

### GEQ-normalised ARG abundance

```
ARG copies per genome equivalent =  RPKM_ARG / GEQ
```

This is the most biologically interpretable ARG abundance metric — it answers: *"For every bacterial genome in this sample, how many copies of this resistance gene are present?"*

---

## 4. Building the KMA database — step by step

See `build_gtdb_scg_kma.sh` for the automated version. The manual steps are:

### Step 1: Extract marker gene archives

```bash
mkdir -p gtdb_scg/reps gtdb_scg/all
cd gtdb_scg/
tar -xzf bac120_marker_genes_reps_r232.tar.gz -C ./
```

### Step 2: Filter marker metadata to high-quality SCGs

```bash
# Keep only markers with ≥85% single-copy AND ≥85% ubiquity
awk -F'\t' 'NR==1 || ($5+0 >= 85 && $6+0 >= 85)' \
  bac120_msa_marker_info_r232.tsv \
  > bac120_hq_markers.tsv

echo "High-quality markers: $(tail -n +2 bac120_hq_markers.tsv | wc -l) / 120"
```

### Step 3: Merge all per-marker FASTA into one file

```bash
# Each genome contributes one sequence per marker
# Header format: >genome_accession
# We normalise this to: >genome_accession|marker_id  for KMA compatibility

mkdir -p gtdb_scg/merged

while IFS=$'\t' read -r marker_id name desc length sc_pct ubiq; do
  [[ "$marker_id" == "Marker Id" ]] && continue
  dir="gtdb_scg/reps/${marker_id}"
  [[ ! -d "$dir" ]] && echo "[WARN] Missing: $dir" && continue

  # Merge all FASTA from this marker, normalising headers
  for fasta in "$dir"/*.fna; do
    [[ ! -f "$fasta" ]] && continue
    # Rewrite headers: >genome~marker_copy → >genome|marker
    awk -v mid="$marker_id" '
      /^>/ {
        split(substr($0,2), a, "~")
        genome = a[1]
        # Skip multi-copy instances (copy index > 1)
        if (a[2] ~ /_[2-9]$/ || a[2] ~ /_[0-9][0-9]+$/) skip=1
        else { skip=0; print ">" genome "|" mid }
        next
      }
      !skip { print }
    ' "$fasta"
  done

done < bac120_hq_markers.tsv \
  >> gtdb_scg/merged/bac120_scg_reps_hq.fna

echo "Total sequences: $(grep -c '^>' gtdb_scg/merged/bac120_scg_reps_hq.fna)"
```

### Step 4: Quality filter — remove multi-copy and short sequences

```bash
# Remove sequences shorter than 75 bp (fragments)
seqkit seq -m 75 gtdb_scg/merged/bac120_scg_reps_hq.fna \
  > gtdb_scg/merged/bac120_scg_reps_hq_len75.fna

echo "After length filter: $(grep -c '^>' gtdb_scg/merged/bac120_scg_reps_hq_len75.fna)"
```

### Step 5: Deduplicate

```bash
seqkit rmdup -j 8 -s -i \
  gtdb_scg/merged/bac120_scg_reps_hq_len75.fna \
  > gtdb_scg/merged/bac120_scg_reps_hq_len75_dedup.fna

echo "After dedup: $(grep -c '^>' gtdb_scg/merged/bac120_scg_reps_hq_len75_dedup.fna)"
```

### Step 6: Build KMA index

```bash
kma index \
  -i gtdb_scg/merged/bac120_scg_reps_hq_len75_dedup.fna \
  -o /export/data/kma_db/bac120_scg_reps_r232 \
  -verbose

# Expected output files:
# bac120_scg_reps_r232.index
# bac120_scg_reps_r232.seq
# bac120_scg_reps_r232.name
# bac120_scg_reps_r232.length.b
# bac120_scg_reps_r232.comp.b
```

**Expected database size:** ~1–3 GB for the index.

---

## 5. KMA alignment strategy

### KMA flags for SCG abundance estimation

```bash
kma \
  -ipe sample_R1.fastq.gz sample_R2.fastq.gz \
  -o sample_scg \
  -t_db /export/data/kma_db/bac120_scg_reps_r232 \
  -t 16 \
  -1t1 \        # one template per read (competitive mapping — critical)
  -mem_mode \   # memory-efficient mode for large databases
  -and \        # both paired reads must map (reduces false positives)
  -apm f \      # paired-end penalty: fragment mode
  -ef \         # extended output with fragment counts
  -sam \        # output SAM for per-read analysis (optional)
  -bcNano       # do NOT use; this is for Nanopore
```

**Critical flag: `-1t1`** — forces each read pair to map to exactly one template (the best match). Without this, KMA distributes reads across all matching templates, which is correct for many use cases but inflates genome-equivalent estimates when multiple species share similar SCG sequences.

### Comparison with ResFinder KMA flags

| Parameter | SCG / GEQ database | ResFinder ARG database | Reason for difference |
|-----------|-------------------|----------------------|----------------------|
| `-1t1` | **Yes** | Yes | Competitive mapping |
| `-mem_mode` | Yes | Yes | Large DB |
| `-and` | Yes | Yes | PE quality |
| `-apm f` | Yes | Yes | PE mode |
| `-ef` | **Yes** | Yes | Extended output needed |
| Identity cutoff | ≥80% | ≥90% | SCGs are more divergent |
| Coverage cutoff | ≥80% | ≥90% | Partial gene hits acceptable |

The lower identity cutoff for SCGs is intentional — you want to capture reads from organisms whose SCGs differ from the representative genome by up to ~20%, which is normal at species-cluster resolution.

---

## 6. Interpreting KMA output for GEQ

### KMA output files

After running `kma`, the output consists of:

| File | Description |
|------|-------------|
| `sample.res` | Per-template summary: identity, coverage, depth, reads |
| `sample.frag.gz` | Per-fragment (read pair) alignment details |
| `sample.mat.gz` | Per-position depth matrix |
| `sample.aln` | Consensus alignment |
| `sample.fsa` | Consensus FASTA |

### The `.res` file — primary input for GEQ

The `.res` file has the same format as ResFinder output:

```
#Template   Score  Expected  Template_length  Template_Identity  Template_Coverage  Query_Identity  Query_Coverage  Depth  q_value  p_value
RS_GCF_000001405.40|PFAM_PF00380.20  1240  1180  121  99.2  100.0  99.2  100.0  12.4  0.0  0.0
RS_GCF_000003025.6|TIGR_TIGR00019   4820  4610  361  97.8  98.9   97.8  98.9   15.1  0.0  0.0
```

The `Template` column encodes `genome_accession|marker_id`, allowing you to:
1. Extract the marker ID → compute per-marker depth
2. Extract the genome accession → link to GTDB taxonomy

### Computing RPKM per SCG marker

```
RPKM_marker = (mapped_reads × 1e6) / (template_length_kb × total_trimmed_reads)
```

In practice, KMA's `Depth` column is already normalised for template length (it is depth per position), so:

```
RPKM_marker ≈ Depth × (1000 / template_length) × (1e6 / total_reads)
```

Or compute directly from `Score` (total mapped bases):

```
RPKM_marker = (Score / template_length) × (1e6 / total_reads)
```

### Computing GEQ

```
GEQ = median(RPKM_marker_i)  for all qualifying markers i
```

The median is preferred over mean because it is robust to occasional markers with anomalously high depth (e.g., from horizontal gene transfer events in a community that happen to share sequences with SCGs).

---

## 7. Per-SCG abundance and taxonomic abundance

### Per-SCG gene-level abundance

This parallels ResFinder's `Template_Identity` approach:
- Each `template_id` in the KMA `.res` file is a specific SCG sequence from a specific genome
- Each `marker_id` (extracted from the template name) groups all templates for one gene

This is analogous to ResFinder where:
- `template_id` = specific allele (e.g. `blaTEM-1_1_AF021808`)
- `gene` = gene family (e.g. `blaTEM`)

For SCGs:
- `template_id` = `RS_GCF_000001405.40|PFAM_PF00380.20` (one genome × one marker)
- `marker` = `PFAM_PF00380.20` (all genomes for that marker)

**Gene-level (marker-level) RPKM** = mean or sum of template-level depths across all genomes for that marker.

**Median across all marker-level RPKMs** = GEQ.

### Taxonomic abundance from SCG hits

The genome accession in each template ID can be joined to `bac120_metadata_r232.tsv.gz` to retrieve GTDB taxonomy. This gives you a taxonomic abundance profile that:

1. Is independent of CCMetagen/Kraken2 (useful for cross-validation)
2. Is anchored to genome-equivalent units (so it's proportional to genome copies, not reads)
3. Uses GTDB taxonomy consistently with the rest of your pipeline

```
Template hit → genome accession → metadata → GTDB taxonomy
RS_GCF_000001405.40 → d__Bacteria;p__Pseudomonadota;...;s__Escherichia coli
```

The depth of each template, normalised by GEQ, gives the relative abundance of each genome (and by extension each taxon) in genome-equivalent units.

---

## 8. Aligning with ResFinder output format

See `combine_gtdb_resfinder.R` for the implementation. The conceptual alignment:

| ResFinder concept | SCG/GTDB equivalent |
|-------------------|---------------------|
| `templateID` | `genome_accession|marker_id` |
| `gene` (family) | `marker_id` (e.g. PFAM_PF00380.20) |
| `RPKM_TemplateID` | RPKM per genome × marker |
| `RPKM_AMRClass` | RPKM per marker (all genomes) |
| `drugClass` | `marker_name` (gene function) |
| `amrClass` | `marker_id` |
| `trimmed.total_reads` | same from sample metadata |
| `GEQ` | median(RPKM per marker) — **new column** |
| `RPKM_ARG / GEQ` | ARGs per genome equivalent — **new column** |

---

## 9. References

| Reference | Relevance |
|-----------|-----------|
| Parks DH, et al. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalised and complete genome-based taxonomy. *Nucleic Acids Research*, 50(D1), D785–D794. https://doi.org/10.1093/nar/gkab776 | GTDB r207+ methods |
| Chaumeil PA, et al. (2022). GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database. *Bioinformatics*, 38(23), 5315–5316. https://doi.org/10.1093/bioinformatics/btac672 | GTDB-Tk tool using bac120 markers |
| Manor O & Borenstein E. (2015). MUSiCC: a marker genes based framework for metagenomic normalization and accurate profiling of gene abundances in the microbiome. *Genome Biology*, 16, 53. https://doi.org/10.1186/s13059-015-0610-8 | GEQ normalisation theory |
| Nayfach S & Pollard KS. (2015). Average genome size estimation improves comparative metagenomics and sheds light on the functional ecology of the human microbiome. *Genome Biology*, 16, 51. https://doi.org/10.1186/s13059-015-0611-7 | Genome size normalisation |
| Clausen PTLC, Aarestrup FM, Lund O. (2018). Rapid and precise alignment of raw reads against redundant databases with KMA. *BMC Bioinformatics*, 19, 307. https://doi.org/10.1186/s12859-018-2336-6 | KMA aligner |
| Bortolaia V, et al. (2020). ResFinder 4.0 for predictions of phenotypes from genotypes. *Journal of Antimicrobial Chemotherapy*, 75(12), 3491–3500. https://doi.org/10.1093/jac/dkaa345 | ResFinder output format reference |
| Hendriksen RS, et al. (2019). Global monitoring of antimicrobial resistance based on metagenomics analyses of urban sewage. *Nature Communications*, 10, 1124. https://doi.org/10.1038/s41467-019-08853-3 | GEQ in wastewater AMR context |

# Pathogen Detection Validation Subworkflow

A Nextflow DSL2 subworkflow for orthogonal validation of Kraken2-detected pathogens via read extraction and BLAST confirmation. Built to nf-core standards.

---

## Table of Contents

- [Purpose](#purpose)
- [Architecture](#architecture)
- [File layout](#file-layout)
- [Inputs and outputs](#inputs-and-outputs)
- [Usage example](#usage-example)
- [Configuration](#configuration)
- [Validation criteria](#validation-criteria)
- [Alternative validation strategies](#alternative-validation-strategies)
- [References](#references)

---

## Purpose

Kraken2 is a k-mer-based classifier optimised for speed. Its taxonomic assignments are accurate at the genus level for most well-represented organisms but can produce false-positive species-level calls due to:

- **K-mer ambiguity** — short k-mers shared between related species
- **Database contamination** — misannotated reference genomes
- **Horizontal gene transfer** — mobile elements producing chimeric taxonomic signals
- **Sample contamination** — laboratory or reagent-introduced reads
- **Low-coverage hits** — single-read assignments to rare taxa

For pathogen surveillance — where false positives can trigger costly public health responses and false negatives can miss real outbreaks — every clinically significant detection should be validated by an orthogonal method. This subworkflow extracts the actual reads Kraken2 assigned to each target pathogen and BLASTs them against a trusted reference database, confirming the assignment with full-length sequence alignment rather than k-mer hashing.

---

## Architecture

```
┌──────────────────────────────────────────────────────────────────────┐
│  INPUT                                                               │
│  - Paired-end FASTQ (clean, host-decontaminated)                     │
│  - Kraken2 classification output (per-read assignments)              │
│  - Kraken2 report file                                               │
│  - List of target pathogen taxids                                    │
│  - Reference BLAST database (local or remote)                        │
└──────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌──────────────────────────────────────────────────────────────────────┐
│  STEP 1: KRAKENTOOLS_EXTRACTKRAKENREADS                              │
│   Tool: KrakenTools/extract_kraken_reads.py (Lu et al. 2017)         │
│   Action: Pulls reads classified to target taxids + descendants      │
│   Output: FASTA (or FASTQ) of extracted reads, one file per sample   │
└──────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌──────────────────────────────────────────────────────────────────────┐
│  STEP 2: CONCAT_EXTRACTED_READS                                      │
│   Action: Combines R1 + R2 FASTA into a single query file            │
│           with /1 /2 suffix on read IDs                              │
└──────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌──────────────────────────────────────────────────────────────────────┐
│  STEP 3: BLAST_BLASTN                                                │
│   Tool: NCBI BLAST+ blastn                                           │
│   Modes: 'local' | 'custom' | 'remote'                               │
│   Output: Tabular hits (outfmt 6) with taxonomy fields               │
│           Per-taxon summary (hits, mean identity, coverage, bitscore)│
└──────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌──────────────────────────────────────────────────────────────────────┐
│  STEP 4: VALIDATE_HITS                                               │
│   Action: Cross-references Kraken2 claim vs BLAST evidence           │
│           Applies thresholds (n_hits, %identity, qcov, bitscore)     │
│   Output: Per-pathogen validation status (CONFIRMED|WEAK|FAILED)     │
└──────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌──────────────────────────────────────────────────────────────────────┐
│  OUTPUTS                                                             │
│  - validation.tsv  : per-sample × per-taxid status table             │
│  - validation.json : same as TSV in machine-readable form            │
│  - blast.tsv       : raw BLAST hits for downstream inspection        │
│  - blast.summary   : aggregated per-taxon BLAST stats                │
└──────────────────────────────────────────────────────────────────────┘
```

---

## File layout

```
your_pipeline/
├── main.nf
├── nextflow.config
├── conf/
│   └── modules.config
├── modules/
│   └── local/
│       ├── krakentools/
│       │   └── extractkrakenreads/
│       │       ├── main.nf
│       │       └── meta.yml
│       └── blast/
│           └── blastn/
│               ├── main.nf
│               └── meta.yml
└── subworkflows/
    └── local/
        └── validate_pathogen_detection/main.nf
```

---

## Inputs and outputs

### Inputs to the subworkflow

| Parameter | Type | Description |
|-----------|------|-------------|
| `ch_input` | channel | tuple `val(meta)`, `path(reads)`, `path(kraken_output)`, `path(kraken_report)` |
| `target_taxids` | list | NCBI taxIDs to extract and validate (e.g. `[562, 1639, 666]`) |
| `blast_db_dir` | path | directory containing BLAST DB (or `[]` for remote mode) |
| `blast_db_name` | string | BLAST DB name (e.g. `'nt'`, `'core_nt'`, `'pathogen_refseq'`) |

### Outputs emitted

| Emit channel | Description |
|--------------|-------------|
| `extracted_reads` | R1 FASTA/FASTQ of extracted reads |
| `extracted_r2` | R2 FASTA/FASTQ (if paired) |
| `extract_summary` | KrakenTools extraction summary table |
| `blast_tsv` | Full BLAST tabular output |
| `blast_summary` | Per-taxon aggregated BLAST stats |
| `validation_report` | Final per-pathogen validation TSV |
| `versions` | nf-core versions.yml |

---

## Usage example

```groovy
// main.nf
nextflow.enable.dsl = 2

include { KRAKEN2 }                       from './modules/local/kraken2/main'
include { VALIDATE_PATHOGEN_DETECTION }   from './subworkflows/local/validate_pathogen_detection'

workflow {

    // 1. Read samplesheet
    ch_reads = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, single_end: false]
            tuple(meta, [file(row.fastq_1), file(row.fastq_2)])
        }

    // 2. Run Kraken2 classification
    KRAKEN2(ch_reads, file(params.kraken2_db))

    // 3. Combine reads with kraken outputs for validation subworkflow
    ch_for_validation = ch_reads
        .join(KRAKEN2.out.classification)
        .join(KRAKEN2.out.report)
        .map { meta, reads, kraken_out, kraken_rep ->
            tuple(meta, reads, kraken_out, kraken_rep)
        }

    // 4. Target pathogen taxids (e.g. Vibrio cholerae, Salmonella Typhi, E. coli)
    target_taxids = [666, 90370, 562, 1639, 1280, 287, 573]

    // 5. Validate
    VALIDATE_PATHOGEN_DETECTION(
        ch_for_validation,
        target_taxids,
        file(params.blast_db_dir),
        params.blast_db_name
    )

    // 6. Aggregate validation reports across samples
    VALIDATE_PATHOGEN_DETECTION.out.validation_report
        .collectFile(
            name:      'all_samples_validation.tsv',
            keepHeader: true,
            storeDir:  "${params.outdir}/validation_summary"
        )
}
```

### Sample `nextflow.config`

```groovy
params {
    samplesheet     = 'samplesheet.csv'
    outdir          = './results'
    kraken2_db      = '/export/data/kraken2/k2_standard_2024'
    blast_db_dir    = '/export/data/ncbi/blast/db/v5'
    blast_db_name   = 'core_nt'
    publish_dir_mode = 'copy'
}

process {

    withName: 'KRAKENTOOLS_EXTRACTKRAKENREADS' {
        ext.fastq_output    = false
        ext.include_children = true
        ext.include_parents  = false
        ext.max_reads        = 100000      // cap per taxid to avoid huge BLAST input
        cpus   = 4
        memory = '8.GB'
        time   = '2.h'
    }

    withName: 'BLAST_BLASTN' {
        ext.db_mode         = 'local'      // 'local' | 'custom' | 'remote'
        ext.task_mode       = 'megablast'
        ext.evalue          = '1e-20'
        ext.perc_identity   = 95
        ext.qcov_hsp        = 80
        ext.max_target_seqs = 5
        cpus   = 16
        memory = '64.GB'
        time   = '24.h'
    }

    withName: 'VALIDATE_HITS' {
        cpus   = 1
        memory = '2.GB'
        time   = '15.m'
    }
}
```

---

## Configuration

### Choosing BLAST database mode

| Mode | Database location | Speed | Use case |
|------|------------------|-------|----------|
| `local` | On HPC filesystem (e.g. `/export/data/bio/ncbi/blast/db/v5/nr`) | Fast | Production surveillance with downloaded DB |
| `custom` | Custom curated FASTA indexed with `makeblastdb` | Fast | Pathogen-focused database (smaller, faster) |
| `remote` | NCBI remote BLAST API | Very slow | One-off lookups, no local DB available |

For routine surveillance of ~1,400 samples, `local` mode with a downloaded `core_nt` database is recommended. Remote mode is impractical at scale due to NCBI rate limits (~3 requests/second without API key).

### Custom pathogen-focused BLAST database

For faster, more specific validation, build a custom BLAST DB containing only your priority pathogen genomes:

```bash
# 1. Download reference genomes for target pathogens via NCBI datasets
datasets download genome accession GCF_000001405.40 GCF_000005845.2 ...
unzip ncbi_dataset.zip

# 2. Concatenate into single FASTA
cat ncbi_dataset/data/*/*.fna > pathogen_refs.fasta

# 3. Index with makeblastdb (preserve taxonomy)
makeblastdb -in pathogen_refs.fasta \
    -dbtype nucl \
    -title 'WBE_priority_pathogens' \
    -parse_seqids \
    -taxid_map taxid_map.txt \
    -out /export/blast/wbe_pathogens
```

Then use it via:

```groovy
ext.db_mode = 'custom'
params.blast_db_dir  = '/export/blast'
params.blast_db_name = 'wbe_pathogens'
```

---

## Validation criteria

The `VALIDATE_HITS` process applies four thresholds (all tunable):

| Threshold | Default | Rationale |
|-----------|---------|-----------|
| `MIN_BLAST_HITS` | 3 | Single reads can match by chance; require ≥3 |
| `MIN_PIDENT` | 95.0% | Strain-level confirmation; lower for divergent pathogens |
| `MIN_QCOVS` | 80.0% | Avoid partial matches that may indicate contamination |
| `MIN_BITSCORE` | 100.0 | Filters short, low-quality alignments |

Status assignment:

- **CONFIRMED**: All 4 thresholds passed → strong orthogonal validation
- **WEAK_EVIDENCE**: 2–3 thresholds passed → flag for manual review
- **FAILED**: 0–1 thresholds passed → likely Kraken2 false positive
- **NOT_VALIDATED**: BLAST returned no hits to the target taxid → likely false positive

---

## Alternative validation strategies

BLAST validation is the most widely used orthogonal check, but for high-confidence pathogen surveillance you should consider stacking multiple methods. Each addresses different failure modes:

### 1. Assembly-based validation (high-confidence)

**Approach:** Extract Kraken2-assigned reads, perform local de novo assembly (SPAdes, MEGAHIT), then BLAST the assembled contigs rather than raw reads.

**Advantages:**
- Longer query sequences (1–10 kb contigs) give far more specific BLAST hits
- Resolves chimeric reads that confuse k-mer classifiers
- Enables confirmation of full gene-of-interest presence (e.g. virulence factors, toxin genes)

**Disadvantages:**
- Requires sufficient read depth (typically ≥30× coverage of the target genome region)
- Slower; assembly can take hours per sample

**Implementation hint:** Add a `SPADES_ASSEMBLY` step between `KRAKENTOOLS_EXTRACTKRAKENREADS` and `BLAST_BLASTN`. Use `--meta` mode for SPAdes when extracted reads represent a mixture.

### 2. Marker gene confirmation (PCR-equivalent)

**Approach:** Map extracted reads against a curated set of pathogen-specific marker genes (toxin genes, virulence genes, species-specific loci).

**Examples:**
- *Vibrio cholerae*: `ctxA`, `ctxB`, `tcpA`, `ompW`
- *Salmonella enterica*: `invA`, `ttr`, `hilA`
- *Mycobacterium tuberculosis*: `IS6110`, `rpoB`
- *Klebsiella pneumoniae*: `khe`, `magA`

**Tool:** Use KMA or `bwa mem` against a marker gene FASTA, then check for full-coverage hits.

**Advantages:**
- Highly specific; equivalent to clinical PCR confirmation
- Fast and computationally cheap
- Identifies pathogen virulotypes, not just presence

**Disadvantages:**
- Requires curation of marker gene database
- May miss novel strains lacking canonical markers

### 3. Bracken abundance re-estimation

**Approach:** Use Bracken (Bayesian Reestimation of Abundance after Classification with Kraken) to re-estimate species-level abundance from Kraken2 reports, addressing Kraken2's tendency to assign reads to higher taxonomic levels when species-level discrimination is ambiguous.

**Tool:** [Bracken](https://github.com/jenniferlu717/Bracken) — same authors as KrakenTools.

**Use case:** Cross-check Kraken2 species calls against Bayesian-corrected abundance estimates. Pathogens detected by Kraken2 but absent from Bracken output (or vice versa) warrant investigation.

### 4. MGS2AMR-style targeted assembly

**Approach:** For pathogens carrying clinically relevant genes (resistance, virulence), use MGS2AMR or similar gene-seeded assembly tools to reconstruct the gene + flanking genomic context, then BLAST the assembled context to confirm host species.

**Advantages:**
- Confirms not just presence of the pathogen, but the gene of interest within that pathogen
- Distinguishes chromosomal vs plasmid-borne ARGs/virulence factors

**Use case:** Already part of your ARG pipeline; can be reused for pathogen validation by seeding assembly with pathogen marker genes instead of ARGs.

### 5. Mapping-based confirmation with strict depth/breadth thresholds

**Approach:** Map all sample reads (not just Kraken2-extracted) against full reference genomes of target pathogens using `bwa mem` or `minimap2`, then require:
- ≥30% breadth of coverage across the genome
- ≥5× mean depth in covered regions
- Even coverage distribution (no extreme localised peaks)

**Tool:** [`covtobed`](https://github.com/telatin/covtobed) or custom scripting on `samtools depth` output.

**Advantages:**
- Validates pathogen presence using the entire genome as evidence, not just k-mers
- Detects very low-abundance pathogens that may be missed by k-mer mappers
- Reveals localised coverage spikes that often indicate contamination or shared mobile elements

**Disadvantages:**
- Slow when run against many reference genomes
- Requires curated reference genome list per pathogen

### 6. Coverage uniformity checks

**Approach:** After mapping extracted reads back to reference, check coverage uniformity. True pathogen presence produces relatively even coverage across the genome; false positives often produce coverage concentrated in conserved/repetitive regions.

**Metric:** Coefficient of variation of per-base depth, or fraction of genome with depth ≥ mean/4.

### 7. Multi-classifier consensus

**Approach:** Run multiple classifiers in parallel (Kraken2, Centrifuge, CCMetagen, MetaPhlAn) and require consensus across ≥2 classifiers before flagging a pathogen as detected.

**Already implemented in your pipeline** — you have Kraken2, Centrifuge, and CCMetagen running. Adding a consensus step that requires agreement across at least two would substantially reduce false positives.

### 8. Negative controls and contamination tracking

**Approach:** Include sample-blank, extraction-blank, and PCR-blank controls in every sequencing run. Subtract or flag any pathogen detected in blanks above background.

**Critical for wastewater surveillance** because reagent contamination with environmental bacteria (especially *Acinetobacter*, *Sphingomonas*, *Burkholderia*) is well documented.

### 9. Strain-level genotyping with StrainGE or sourmash

**Approach:** For confirmed pathogens, use [StrainGE](https://github.com/broadinstitute/StrainGE) or [sourmash](https://github.com/sourmash-bio/sourmash) to assign strain-level identity, enabling outbreak strain tracking across sites and time points.

**Use case:** If *V. cholerae* is detected in multiple Kisumu sites within the same week, is it the same strain (suggesting a point source outbreak) or different strains (suggesting endemic circulation)?

### 10. Recommended validation stack for production WBE

For your wastewater AMR/pathogen surveillance, I recommend layering three methods:

```
1. KRAKENTOOLS extract + BLASTN          ← this subworkflow (fast, broad)
2. SPAdes assembly + BLAST contigs        ← for any CONFIRMED pathogen with ≥100 reads
3. Marker gene mapping (KMA + curated DB) ← virulotype confirmation for top-priority pathogens
```

Only pathogens passing all three layers should be reported as confirmed in surveillance reports. Pathogens passing only the first layer should be flagged for follow-up sequencing or culture confirmation.

---

## References

| Reference | Tool / Method |
|-----------|---------------|
| Lu J, Breitwieser FP, Thielen P, Salzberg SL. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104. https://doi.org/10.7717/peerj-cs.104 | Bracken / KrakenTools |
| Wood DE, Lu J, Langmead B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257. https://doi.org/10.1186/s13059-019-1891-0 | Kraken2 |
| Camacho C, Coulouris G, Avagyan V, et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. https://doi.org/10.1186/1471-2105-10-421 | BLAST+ |
| Marcelino VR, Clausen PTLC, Buchmann JP, et al. (2020). CCMetagen: comprehensive and accurate identification of eukaryotes and prokaryotes in metagenomic data. *Genome Biology*, 21, 103. https://doi.org/10.1186/s13059-020-02014-2 | CCMetagen |
| Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. (2017). metaSPAdes: a new versatile metagenomic assembler. *Genome Research*, 27, 824–834. https://doi.org/10.1101/gr.213959.116 | metaSPAdes |
| van Dijk LR, Walker BJ, Straub TJ, et al. (2022). StrainGE: a toolkit to track and characterize low-abundance strains in complex microbial communities. *Genome Biology*, 23, 74. https://doi.org/10.1186/s13059-022-02630-0 | StrainGE |
| Pierce NT, Irber L, Reiter T, Brooks P, Brown CT. (2019). Large-scale sequence comparisons with sourmash. *F1000Research*, 8, 1006. https://doi.org/10.12688/f1000research.19675.1 | sourmash |
| Ewels PA, Peltzer A, Fillinger S, et al. (2020). The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology*, 38, 276–278. https://doi.org/10.1038/s41587-020-0439-x | nf-core standards |

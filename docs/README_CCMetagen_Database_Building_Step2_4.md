## Why 300 is the correct minimum for your data

### The case for MIN_LEN=300 with 150×150 paired-end reads

With 150×150 PE reads, each read is 150 bp. KMA processes paired-end reads either as individual reads or as a fragment (the insert spanning both reads). The minimum template length should satisfy **both** alignment geometry and statistical discriminability.

**Argument 1 — Template must be longer than a single read**

When a reference sequence is shorter than the read being mapped, the read overhangs the template ends. KMA handles this with end-gap penalties, but partial alignments covering only a fraction of the reference are unreliable. A 150 bp reference aligned by a 150 bp read provides no flanking context — the aligner cannot distinguish whether the read belongs to this template or to one of dozens of similar short sequences in the database.

At MIN_LEN=150, every template is exactly one read long. At MIN_LEN=300, every template is at least two read lengths long, meaning even a single read can be fully contained within the template with flanking sequence on at least one side.

**Argument 2 — Paired-end fragment coverage requires template ≥ insert size**

In paired-end sequencing at 150×150, the typical insert size is 300–500 bp. When KMA processes a paired-end fragment in `-apm f` (fragment) mode, it expects to place both reads on the same template simultaneously. If the template is shorter than the insert size, at least one read will overhang or not place at all. A 300 bp minimum ensures the template is long enough to accommodate at least one complete read pair at the shortest realistic insert size.

**Argument 3 — Statistical discrimination requires sufficient unique k-mer content**

KMA builds a k-mer index (default k=16 for nucleotide databases). A 150 bp sequence contains 135 unique 16-mers. A 300 bp sequence contains 285 unique 16-mers. Doubling the sequence length nearly doubles the number of discriminatory k-mers available to distinguish this template from others. Shorter sequences are more likely to share all their k-mers with other templates, leading to multi-mapping and abundance inflation.

For the `nt` and `core_nt` databases, which contain millions of sequences including highly similar variants differing by 1–2 SNPs, this k-mer discriminability argument is particularly strong. A 150 bp sequence might share all 135 of its k-mers with a near-identical variant; a 300 bp sequence is far less likely to have complete k-mer overlap with any other template.

**Argument 4 — 200 bp is too close to read length to add meaningful benefit**

The case for 200 over 300 is essentially: "we want to keep more sequences." But the sequences lost between 200 and 300 bp are short fragments that are:
- Too short to be reliably classified by paired-end mapping
- Disproportionately likely to be partial gene submissions, PCR products, or sequencing artefacts
- Absent from complete reference genomes (which are the true targets of CCMetagen classification)

The gain from keeping 200–300 bp sequences is marginal recovery of a few fragmented or historical submissions. The cost is increased multi-mapping noise, slower KMA indexing, and a larger database with more ambiguous k-mer entries.

**Concrete summary:**

| Threshold | Rationale | Risk |
|-----------|-----------|------|
| 150 bp | Matches read length | Every template fully overhangs; no discriminating flanking context |
| 200 bp | Slightly longer than read | Still too short for PE fragment placement; minimal k-mer gain |
| **300 bp** | **Two read lengths; minimum PE insert** | **Optimal balance of sensitivity and specificity** |
| 500 bp | Conservative | Loses some legitimate short viral/plasmid sequences |

---

## How the dustmasker code works

Breaking the pipeline down step by step:

```bash
dustmasker -in "$STEP4_INPUT" -outfmt fasta \
  | awk '
      /^>/ {print; next}
      {gsub(/[a-z]/,"N"); print}
    ' \
  > "$STEP4"
```

### Stage 1: `dustmasker -in input.fasta -outfmt fasta`

`dustmasker` implements the **DUST algorithm** (Morgulis et al. 2006), which identifies low-complexity regions by computing a score based on triplet frequencies within a sliding window.

For each position in the sequence, it calculates:

```
score = (sum of triplet_count × (triplet_count - 1) / 2) / (window_length - 2)
```

Where `triplet_count` is the count of each distinct 3-mer in the window. A window dominated by a single repeated triplet (e.g. `ATATAT...`) has a very high score; a window with even triplet distribution has a low score. The default threshold is 2.5 — windows scoring above this are considered low-complexity.

The output is FASTA where low-complexity regions are **softmasked** — converted to lowercase letters while the rest remains uppercase:

```
>9606|Homo sapiens chromosome 1
ATCGATCGatatatatatATATATATGCGCGCGCGcgcgcgcgATCGATCG
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         softmasked (lowercase) = low-complexity
```

High-complexity regions remain uppercase (`ATCGATCG`). The sequence content is preserved — no information is lost at this stage.

### Stage 2: The `awk` hard-masking conversion

```awk
/^>/ {print; next}         # print header lines unchanged
{gsub(/[a-z]/,"N"); print} # replace every lowercase letter with N
```

This converts the softmasked output to **hard masking** — lowercase letters become `N` (unknown nucleotide):

```
Before: ATCGATCGatatatatatATATATATGCGCGCGCGcgcgcgcgATCGATCG
After:  ATCGATCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATCGATCG
```

Hard masking is irreversible. The masked regions are completely invisible to any downstream aligner — KMA cannot map reads to `N` positions, cannot build k-mers from them, and effectively ignores them.

---

## Is dustmasking necessary for CCMetagen?

**Short answer: No, and for your specific workflow it is recommended to set `MASK_LOW_COMPLEX=0`.**

Here is the detailed reasoning:

### What CCMetagen/KMA already does about low-complexity

KMA's k-mer indexing naturally handles repetitive sequences through two mechanisms:

**1. Multi-mapping k-mers receive lower weight.** K-mers that appear in many templates are down-weighted in the alignment score because they carry less discriminatory information. A k-mer from a poly-A tail that appears in 50,000 templates contributes almost nothing to template selection. This is functionally equivalent to soft-masking — repetitive k-mers are deprioritised without being discarded.

**2. The `-1t1` competitive mapping flag.** This ensures each read maps to exactly one template. Even if a read partially matches a low-complexity region shared by many templates, KMA picks the one with the best overall alignment score across the full template. The competitive nature of this mapping means low-complexity matches are penalised relative to specific matches.

### The three risks of hard masking for a classification database

**Risk 1 — You destroy diagnostic signal in primers and conserved regions**

Many taxonomically diagnostic sequences contain short low-complexity stretches flanked by highly specific regions. For example, bacterial 16S variable regions contain homopolymers and short repeats that DUST would mask, but these are real biological sequence contributing to species discrimination. Hard masking converts these to `N`, preventing any read from matching across the masked boundary.

**Risk 2 — Viral genomes are disproportionately affected**

RNA viruses, particularly +ssRNA viruses like Hepatitis C (which you have in your `core_nt` export), have AU-rich 3' regions, internal ribosome entry sites with structured repeats, and other low-complexity features that are genuinely taxonomically informative. Masking these destroys viral classification sensitivity — precisely the pathogens you are most interested in detecting in your wastewater surveillance.

**Risk 3 — Short templates become entirely N after masking**

A 300 bp template that happens to contain a 250 bp low-complexity region (not uncommon for transposable element fragments, satellite sequences, or microsatellite-flanking genomic clones) becomes mostly `N`. KMA cannot build any useful k-mers from it. The sequence survives the length filter but contributes nothing to classification — occupying database space while providing no information. It would have been better to filter it out entirely.

### When masking is genuinely useful

Dustmasking is beneficial in a narrower context:

- **BLAST searches** where ungapped seed extension through low-complexity regions causes spurious high-scoring alignments (`-dust yes` in BLAST+)
- **De novo assembly** databases where repetitive sequences cause misassembly
- **Alignment-based SNP calling** where repetitive regions produce mapping artefacts

KMA's k-mer approach avoids most of these problems by design. The masking step was appropriate for BLAST-era database construction and is now carried forward by convention rather than necessity.

### Recommendation for your pipeline

```bash
# For CCMetagen / KMA databases:
MASK_LOW_COMPLEX=0   # skip dustmasker entirely

# If you want to remove genuinely low-complexity sequences entirely
# (rather than masking within them), use seqkit's complexity filter instead:
seqkit fx2tab -l -g -n input.fasta \
  | awk '$4 > 0.3' \       # keep sequences with GC content > 30% (removes poly-A tracts)
  | cut -f1 \
  | seqkit grep -f - input.fasta \
  > output_complexity_filtered.fasta
```

Or use `prinseq` which calculates DUST score per sequence and lets you filter whole sequences below a threshold rather than masking within sequences:

```bash
prinseq-lite.pl \
  -fasta input.fasta \
  -out_good stdout \
  -out_bad /dev/null \
  -lc_method dust \
  -lc_threshold 7    # remove sequences where >7% of sequence is low-complexity
```

This filters out genuinely low-quality sequences (poly-homopolymer submissions, vector-only sequences that passed title filtering) without destroying diagnostic signal within legitimate genomic sequences.

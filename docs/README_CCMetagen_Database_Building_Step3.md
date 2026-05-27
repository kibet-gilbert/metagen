# How does deduplication affect the `ccmetagen` KMA database?

This is an important question that touches on a real tension in database construction. Let me work through it carefully.

---

## What `seqkit rmdup -s -i` actually does

The `-s` flag deduplicates by **exact sequence identity** — two sequences must be byte-for-byte identical to be collapsed. The `-i` flag makes the comparison case-insensitive (treating `A` and `a` as the same). One representative is kept; all others are removed regardless of their headers, taxids, or organism of origin.

The stated rationale — "identical sequences from different submitters inflate database size" — is correct but incomplete. The real picture is more nuanced.

---

## Does deduplication affect classification? Yes, significantly

### The core problem: KMA is a competitive mapper

KMA assigns each read to the **single best-matching template** in the database (`-1t1` mode). When multiple identical sequences exist in the database with different headers (different taxids), KMA must choose one arbitrarily or distribute the score. After deduplication, only one representative remains. This has two consequences:

**Consequence 1 — The surviving taxid may be wrong for your sample**

Suppose the sequence `ATCGATCG...` (300 bp, a 16S hypervariable region) exists in the database as:
- `>1280|Staphylococcus aureus strain X`
- `>1282|Staphylococcus epidermidis strain Y`
- `>562|Escherichia coli strain Z`

After deduplication, only one survives — say `1280`. Every read matching this sequence in your wastewater sample now gets attributed to *S. aureus*, even if the true source is *E. coli* or *S. epidermidis*. The classification is not wrong in the sense that the sequence matches — it is wrong in the sense that the taxid assignment is arbitrary.

This matters enormously for highly conserved loci. The GTDB bac120 markers, certain 16S regions, many housekeeping genes, and insertion sequences are genuinely shared across very divergent taxa. Deduplication collapses real biological diversity into a single arbitrary representative.

**Consequence 2 — Abundance estimates are not affected by deduplication itself**

Here is where the concern is partially addressed: if your goal is **quantification** (RPKM, GEQ) rather than **species identification**, deduplication of identical sequences does not change the total read count. Reads that would have mapped to three identical sequences now map to one. The depth on that one template is the same as the combined depth would have been, because KMA competitive mapping would have assigned each read to exactly one template anyway.

So for GEQ normalisation specifically — where you care about read depth per SCG marker, not which species the marker came from — deduplication is relatively safe.

---

## The concerns in full

### 1. Highly conserved sequences lose taxonomic resolution (partially unaddressed)

You identify this correctly. A sequence present in 500 genomes across 50 genera gets collapsed to one. The surviving header is essentially random (seqkit keeps the first occurrence in file order). For `nt` and `core_nt`, the file order depends on how blastdbcmd iterates the database — not on biological priority.

**What is not addressed:** the script makes no attempt to choose a *representative* header based on taxonomy. It keeps an arbitrary one. For a CCMetagen database used for taxonomic classification, this means:
- Reads from organisms whose sequence was collapsed will be attributed to whichever organism happened to survive
- You lose all information about the true distribution of that sequence across taxa

**Partial mitigation:** KMA's output includes `Template_Identity` and `Template_Coverage`. Very conserved sequences with many cross-taxon hits will tend to show low specificity in the KMA results anyway, and CCMetagen applies its own post-processing to handle ambiguous hits. But the deduplication step removes information that CCMetagen could otherwise use to resolve ambiguity.

### 2. Sequences that are identical at nucleotide level but encode different taxa are not the same biologically (unaddressed)

Two sequences can be 100% identical at nucleotide level but represent genuinely different lineages — for example, due to horizontal gene transfer, convergent evolution, or insufficient resolution in a short marker. Deduplication treats these as redundant when they are not.

**What is not addressed:** the script has no taxonomy-aware deduplication. A better approach for classification databases is to cluster at 95–99% identity using `vsearch` or `usearch` while keeping the full cluster membership in a lookup table, so that a hit to the centroid can be expanded back to all member taxa for LCA (lowest common ancestor) resolution — exactly what tools like Kraken2 and DIAMOND do internally.

### 3. Deduplication by exact sequence identity misses near-identical redundancy (partially unaddressed)

Two sequences that differ by one SNP are not deduplicated. For `nt`, thousands of sequences differ by 1–2 bases due to sequencing error in the original submissions. These are biologically meaningless variants that still inflate the database.

**What is addressed:** the length filter (≥300 bp) partially mitigates this by removing very short fragments that are most likely to be spurious near-identical variants.

**What is not addressed:** a clustering step at 99% identity would be more thorough. This is the approach used by SILVA for 16S databases and by UniRef for protein databases.

### 4. The duplicate list is written but not used (unaddressed)

The `duplicates.txt` file records what was removed, but there is no downstream step that uses this information. If a collapsed sequence belongs to a pathogen you care about (e.g., a conserved region of *Vibrio cholerae* that is also found in environmental *Vibrios*), you have no way of knowing from the current pipeline that the pathogen-specific entry was the one removed.

---

## What would fully address these concerns

For a CCMetagen classification database, the correct approach is **taxonomy-aware clustering** rather than exact deduplication:

```bash
# Step 1: Cluster at 95% identity, keeping all members
vsearch \
  --cluster_fast input.fasta \
  --id 0.95 \
  --centroids centroids.fasta \
  --uc clusters.uc \
  --threads 16

# Step 2: Build a cluster membership table: centroid → all member taxids
awk '$1=="H" || $1=="S" {
  centroid = ($1=="S") ? $9 : $10
  member   = $9
  print centroid "\t" member
}' clusters.uc > cluster_membership.tsv

# Step 3: For each centroid, assign LCA taxonomy from member taxids
# Then use this LCA-assigned FASTA as the KMA reference
```

This is how Kraken2 builds its database — sequences sharing a node in the taxonomy tree are assigned to their LCA, so reads mapping to a conserved sequence get attributed to the correct higher-level taxon rather than an arbitrary species. CCMetagen performs its own LCA step post-KMA, but it can only work with the taxids that survive in the database.

---

## Summary

| Concern | Addressed by current script? | Better approach |
|---------|------------------------------|-----------------|
| Database size reduction | ✓ Yes | — |
| Arbitrary taxid survival | ✗ No | Taxonomy-aware LCA assignment |
| Near-identical (99%) redundancy | ✗ Partially | vsearch clustering at 95–99% |
| Loss of pathogen-specific sequences | ✗ No | Priority-list before dedup |
| Audit trail of what was removed | ✓ Partially (duplicates.txt) | Full cluster membership table |
| GEQ quantification accuracy | ✓ Yes (unaffected) | — |
| Classification accuracy for conserved loci | ✗ No | LCA-aware database construction |

For your use case — where CCMetagen is used for community profiling and the GTDB SCG database is used for GEQ quantification — the exact deduplication step is acceptable for the GEQ workflow but suboptimal for the CCMetagen classification workflow. The GEQ analysis does not depend on which taxid a marker is attributed to, only on the total depth across markers. The CCMetagen database would benefit from taxonomy-aware clustering with LCA assignment, which is a substantially more complex build process.

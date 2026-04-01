# **Buildind CCMetagen database from *core_nt* or*nt* NCBI BLAST database**

---

## **Context:**  
We need a ccmetagen database built from NCBI's nt or core_nt databases. It will be used in analysing sequence data containing both DNA and RNA from all organisms. While preparing the database we need to exclude redundant sequences

 - Redundancy removal requires **sequence‑level or taxonomy‑level filtering**
 - With DNA + RNA from all organisms, your combined nt/core\_nt export will include:

*   Genomic DNA
*   Transcript variants (RNA)
*   Predicted mRNAs (XM, XR, XP)
*   Alternative isoforms
*   Pseudogenes
*   Partial fragments
*   Redundant clones
*   Overlapping assemblies from different submitters
*   Contaminants

These **explode redundancy**, slow KMA, and inflate false positives.

---

## **Redundancy should be removed using:**

### ✔ **Sequence-level filtering**

*   **deduplication** (exact): `seqkit rmdup -s -i`
*   **clustering** (90–100%): `vsearch --cluster_fast`
*   **low-complexity masking**: `dustmasker`, `sdust`

### ✔ **Taxonomy-level filtering**

*   Exclude host clades (e.g., Vertebrata, Mammalia)
*   Exclude predicted mRNAs (`PREDICTED:`)
*   Exclude unclassified/environmental sequences (`uncultured`, `metagenomic`, etc.)

---

## **What sequences should you exclude when building a CCMetagen database?**

This is crucial and should exclude:

### **A. Host sequences** (unless you want host hits)

Examples you gave:

    >PREDICTED: Peromyscus californicus...
    >Giant Panda satellite DNA
    >Bos taurus mRNA for bone Gla protein

These are **non‑microbial**, cause **false positives**, and waste disk/compute.

### **B. “PREDICTED:” entries (XM\_, XR\_, XP\_)**

These are computational mRNA/protein predictions:

*   Redundant isoforms
*   Partial transcripts
*   Noisy for mapping
*   Inflate multi-mappers

### **C. Environmental / unclassified / “uncultured” sequences**

*   `uncultured bacterium`
*   `metagenome`
*   `environmental sample`

These add noise and non‑interpretable hits.

### **D. Partial fragments & tiny sequences**

Filter out sequences <200–300 nt.

### **E. Patents, vectors, adaptors**

Often labeled as:

*   “vector”
*   “clone”
*   “synthetic construct”
*   “patented DNA”

---

## **Recommended workflow for a high-quality CCMetagen database**

### **Step 1: Export from nt/core\_nt**

(Use fixed outfmt string)

```bash
blastdbcmd -db nt \
  -entry all \
  -outfmt $">%T|%t\n%s" \
  | pigz > nt.raw.fasta.gz
```

### **Step 2 (Optional): Taxonomy-based filtering**

Exclude:

*   Vertebrata (7742)
*   Metazoa (33208)
*   Embryophyta (3193)
*   Any host species

Using taxonkit or direct lists:

```bash
blastdbcmd -db nt \
  -entry all \
  -negative_taxidlist host_exclude_list.txt \
  -outfmt $">%T|%t\n%s" \
  > nt.microbes.fasta
```

### **Step 3: Title filtering**

Remove predicted entries:

```bash
grep -v -i '^>[^|]*|PREDICTED:' nt.microbes.fasta > nt.clean_titles.fasta
```

Remove environmental:

```bash
grep -v -Ei 'uncultured|environmental|metagenome' nt.clean_titles.fasta > nt.clean2.fasta
```

Remove vectors/synthetic constructs:

```bash
grep -v -Ei 'vector|synthetic|cloning' nt.clean2.fasta > nt.clean3.fasta
```

### **Step 4: Length filter**

```bash
seqkit seq -m 200 nt.clean3.fasta > nt.len200.fasta
```

### **Step 5: Deduplicate**

```bash
seqkit rmdup -s -i nt.len200.fasta > nt.dedup.fasta
```

### **Step 6: Compress & index for KMA**

```bash
kma index -i nt.dedup.fasta -o nt.kma
```

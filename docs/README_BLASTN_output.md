## Column-by-column breakdown

Using your first hit as the example:

```
qseqid   = 1/VH01928:52:222CHTKNX:1:1101:57678:1120/1
sseqid   = gi|2034578690|gb|CP048740.1|
pident   = 97.987
length   = 149
mismatch = 3
gapopen  = 0
qstart   = 1
qend     = 149
sstart   = 2380846
send     = 2380994
evalue   = 1.06e-64
bitscore = 259
staxids  = 2708052
sscinames= N/A
scommons = N/A
qcovs    = 99
qcovhsp  = 99
stitle   = Romboutsia sp. 13154 chromosome
```

---

### Query identity fields

**`qseqid`** — your read ID. The `1/` prefix is from the combined R1+R2 FASTA you built in `CONCAT_EXTRACTED_READS` (the `i=1` counter prepended to each read header). `VH01928:52:222CHTKNX:1:1101:57678:1120` is the Illumina read name. `/1` means it's the R1 mate.

**`qstart` / `qend`** — which bases of YOUR read were aligned. Here `1–149` out of a 150bp read, so essentially the full read was used.

**`qcovs`** — query coverage across all HSPs (High-Scoring Pairs) for this subject. `99` means 99% of your read aligned to this reference sequence. Above 80% is generally considered a good alignment.

**`qcovhsp`** — query coverage for this single HSP only. Here it equals `qcovs` because there's only one alignment block, no gaps or fragmentation.

---

### Subject (reference) identity fields

**`sseqid`** — the matched reference sequence in the database. `gi|2034578690|gb|CP048740.1|` breaks down as:
- `gi|2034578690|` — legacy NCBI GI number
- `gb|CP048740.1|` — GenBank accession CP048740, version 1

**`sstart` / `send`** — coordinates on the reference where your read aligned. `2380846–2380994` means your read hit at position 2.38 million on a chromosome. When `sstart > send` (as in hit 3: `15572–15423`), the alignment is on the **reverse strand** — the read matched the reverse complement.

**`staxids`** — NCBI taxonomy ID of the hit organism. `2708052` = *Romboutsia sp. 13154*. You can look this up at `https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2708052`. Multiple taxids would appear semicolon-separated if a sequence maps to multiple taxa.

**`sscinames` / `scomnames`** — scientific and common names. Both show `N/A` here because your BLAST database was built without the taxonomy name files (`taxdb.btd` / `taxdb.bti`). You're getting taxids correctly but not the string names. Fix:
```bash
# Download taxonomy name files into your blast DB directory
cd /export/data/bio/ncbi/blast/db/v5
update_blastdb.pl taxdb
# or wget directly:
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -xzf taxdb.tar.gz
```

**`stitle`** — free-text description from the sequence record. `Romboutsia sp. 13154 chromosome` tells you what the matched sequence is.

---

### Alignment quality fields

**`pident`** — percentage of aligned positions that are identical matches. `97.987%` over 149bp means 3 mismatches, which is excellent. Interpretation:

| pident | Interpretation |
|---|---|
| ≥ 99% | Near-identical — likely same strain |
| 97–99% | Same species, minor variation |
| 95–97% | Same species, strain divergence |
| 90–95% | Same genus, different species |
| < 90% | Family level or below — treat with caution |

**`length`** — number of base pairs in the alignment (including mismatches and gaps, not just identical bases). `149` out of a 150bp read = almost the full read aligned.

**`mismatch`** — number of non-identical aligned positions. `3` mismatches over 149bp = the 2.013% non-identity in pident.

**`gapopen`** — number of gap openings. `0` = no insertions or deletions, a clean ungapped alignment.

**`evalue`** — probability that this hit is a random match. `1.06e-64` means there is a 1 in 10⁶⁴ chance this match is random — essentially impossible. Rule of thumb:

| evalue | Significance |
|---|---|
| < 1e-50 | Highly significant |
| 1e-50 to 1e-10 | Significant |
| 1e-10 to 0.01 | Borderline — check manually |
| > 0.01 | Not significant |

**`bitscore`** — normalised alignment score independent of database size. Higher is better. `259` is strong. Use this to rank hits when evalues are similarly tiny (e.g. comparing `1.06e-64` vs `1.37e-63`).

---

## Reading the four hits together

```
Hit 1: pident=97.99%, len=149, evalue=1.06e-64, bitscore=259 → Romboutsia sp. 13154
Hit 2: pident=95.97%, len=149, evalue=1.07e-59, bitscore=243 → Romboutsia sp. 13154  (same accession, different position on chromosome)
Hit 3: pident=97.33%, len=150, evalue=1.37e-63, bitscore=255 → Paenalkalicoccus suaedae M4U3P1  (reverse strand — sstart > send)
Hit 4: pident=97.33%, len=150, evalue=1.37e-63, bitscore=255 → Paenalkalicoccus suaedae  (different position, same score as hit 3)
```

**What this means biologically:** this read aligns well to two different organisms (*Romboutsia sp.* and *Paenalkalicoccus suaedae*). This is expected for a metagenomics read — it either:
- Comes from a conserved region (e.g. 16S rRNA, housekeeping gene) present in multiple organisms, or
- Is genuinely ambiguous and the true source organism can't be determined from this read alone

`validate_pathogen_detection` subworkflow downstream `VALIDATE_HITS` step handles exactly this. It cross-references Kraken2's classification with these BLAST hits and counts how many reads support each taxid, filtering out ambiguous assignments using `min_blast_hits` and `min_pident` thresholds.

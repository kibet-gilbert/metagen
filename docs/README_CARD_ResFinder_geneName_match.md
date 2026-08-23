# CARD &harr; ResFinder gene-name reconciliation: methodology, validation, and citations

**Task**: build a merged dictionary matching CARD's Antibiotic Resistance Ontology (ARO)
terms (from `card_prevalence.txt`) to ResFinder DB acquired-gene names (from `notes.txt`,
enriched with `phenotypes.txt`), for downstream comparative prevalence analysis.

**Result**: 1,026 of 2,713 CARD ARO terms (37.8%) matched to a ResFinder gene name with
documented, rule-based provenance. Every match was produced by a deterministic,
re-runnable rule, not free-form guessing; the full pipeline is in `build_card_resfinder_mapping.py`.

---

## 0. Reproducing this analysis (CLI usage)

Both scripts are self-documenting (`--help`); this section is a quick reference.

### `build_card_resfinder_mapping.py` -- builds the merged dictionary

```bash
# Defaults: reads /mnt/user-data/uploads, writes /mnt/user-data/outputs
python3 build_card_resfinder_mapping.py

# Custom directories
python3 build_card_resfinder_mapping.py --inputdir ./data --outputdir ./results

# Individually renamed/relocated input files
python3 build_card_resfinder_mapping.py \
    --card-prevalence-file /data/card_v3.9/prevalence.tsv \
    --resfinder-notes-file /data/resfinder_db/notes.txt \
    --resfinder-phenotypes-file /data/resfinder_db/phenotypes.txt \
    --outputdir ./results

# Audit what the safety-net rules actually change
python3 build_card_resfinder_mapping.py --no-category-check --no-manual-overrides -v
```

Full option reference: `python3 build_card_resfinder_mapping.py --help`. Notable flags: `--no-category-check`
disables the drug-class consistency re-ranking (see "Cross-class false positive" in
Section 3 below); `--no-manual-overrides` disables the one literature-curated tie-break;
`--save-intermediate` also writes pickled intermediate DataFrames under
`<outputdir>/intermediate/`; `-v/--verbose` prints every individual override/flag
decision, `-q/--quiet` suppresses non-error output.

### `summarize_mapping.R` -- generates the summary plots below

```bash
# Defaults: reads card_resfinder_merged_dictionary.csv from /mnt/user-data/outputs,
# writes PNGs to /mnt/user-data/outputs/plots
Rscript summarize_mapping.R

# Custom paths, PDF output, larger figures
Rscript summarize_mapping.R --input results/merged.csv \
    --outdir results/plots --format pdf --width 10 --height 6
```

Requires R plus the optparse/readr/dplyr/tidyr/stringr/forcats/ggplot2/scales/patchwork
packages; the script attempts `install.packages()` for anything missing. On a machine
without CRAN access, install the Debian/Ubuntu binaries instead:
`sudo apt-get install r-cran-ggplot2 r-cran-dplyr r-cran-tidyr r-cran-readr r-cran-stringr r-cran-forcats r-cran-scales r-cran-optparse r-cran-patchwork`.

Full option reference: `Rscript summarize_mapping.R --help`.

## 0b. What the plots show

| File | Question it answers |
|---|---|
| `01_match_status_overview` | How many of CARD's 2,713 ARO terms matched, vs. genuinely unmatched, vs. categorically out of ResFinder's acquired-gene scope? |
| `02_match_tier_breakdown` | Which rule (exact / bla-prefix / van-cluster / alias-bridge / ...) found each of the 1,026 matches, and at what confidence? |
| `03_resfinder_class_distribution` | Which ResFinder resistance classes account for the overlap (spoiler: Beta-Lactamase, by a wide margin)? |
| `04_drugclass_coverage` | For CARD's most common drug/mechanism categories, what fraction has a ResFinder counterpart? |
| `05_prevalence_by_match_status` | Does a gene's real-world NCBI prevalence relate to whether it's in ResFinder's scope? (It does -- see below.) |
| `06_summary_dashboard` | Panels 1/2/3/5 combined into one at-a-glance figure. |

**Notable pattern surfaced by panel 5**: the three "out of scope" categories
(chromosomal point mutations, rRNA mutations, regulatory/overexpression) sit at
*substantially higher* NCBI WGS prevalence than either matched or unmatched acquired
genes. This makes biological sense and is worth stating explicitly: a chromosomal
point-mutation "gene" (e.g. *gyrA*) is a variant of a core gene present in essentially
every isolate of a species, so RGI's detection rate approaches 100% by construction,
whereas acquired genes are only present in the subset of strains that received them by
horizontal transfer. It's a reminder to keep those two prevalence notions distinct in
the next analysis phase, rather than pooling all "prevalence" numbers together.

## 1. Why the two databases don't line up 1:1 (scope)

CARD's Resistance Gene Identifier (RGI) annotates five model types: protein homolog
(horizontally acquired genes), protein variant (chromosomal point mutations, e.g. *gyrA*
fluoroquinolone resistance), rRNA gene variant, protein overexpression (regulatory
mutations, e.g. *marR*, *soxS*), and gene cluster models. This reflects CARD's design as
a unified ontology over the whole resistome (McArthur et al. 2013; Jia et al. 2017;
Alcock et al. 2020, 2023).

ResFinder, by contrast, is deliberately split into two companion tools: **ResFinder**
covers only horizontally acquired resistance genes (Zankari et al. 2012; Bortolaia et al.
2020), while chromosomal point mutations are the separate **PointFinder** database
(Zankari et al. 2017) — not provided in this upload. The four files supplied
(`notes.txt`, `phenotypes.txt`, `antibiotic_classes.txt`, `phenotype_panels.txt`) are all
from the acquired-gene side.

**Empirical confirmation from the data itself**: of `card_prevalence.txt`'s 2,713 unique
ARO terms, 2,573 are `protein homolog model` (candidates for a ResFinder match) and 140
are `protein variant model` / `rRNA gene variant model` / `protein overexpression model`
(chromosomal/regulatory — categorically outside ResFinder's scope). Running the full
matching pipeline produced **zero** matches against any of those 140 non-homolog entries
— a clean sanity check that the pipeline isn't cross-contaminating categories.

`phenotype_panels.txt` further confirms ResFinder's species panels are *Campylobacter*,
*E. coli*, *Enterococcus faecalis/faecium*, *Salmonella*, *S. aureus*, and
*M. tuberculosis* — notably **not** *Pseudomonas aeruginosa* or *Acinetobacter
baumannii*. This directly explains the largest single block of CARD-only entries: 491 of
the 1,547 unmatched `protein homolog model` terms are numbered alleles of intrinsic
chromosomal cephalosporinases (PDC = *Pseudomonas*-derived cephalosporinase, ADC =
*Acinetobacter*-derived cephalosporinase, OXY/LEN = *Klebsiella* intrinsic
beta-lactamases) — organisms ResFinder's panels don't cover, so CARD's much finer
allele-level curation (hundreds of PDC/ADC variants) has no ResFinder counterpart by
design, not by matching failure.

## 2. Matching rules, in priority order

All rules were derived empirically from the actual gene-name lists in the two files
(not assumed), then cross-checked against the literature below.

| Tier | Rule | n | Example |
|---|---|---|---|
| `exact` | identical string | 133 | `aadA22` = `aadA22` |
| `case_insensitive` | case-fold only | 164 | `QnrB3` = `qnrB3` |
| `van_structural` | CARD `"vanX gene in vanY cluster"` &harr; ResFinder `"VanX-Y"` | 38 | `vanS gene in vanB cluster` = `VanS-B` |
| `loose` | strip spaces/parens/hyphens, case-fold (apostrophes **preserved**) | 37 | `mphC` = `mph(C)` |
| `debla` | strip ResFinder's leading `bla` prefix | 651 | `CMY-97` = `blaCMY-97` |
| `alias_bridge` | linked only via a `phenotypes.txt` "Alternative name" cross-reference | 26 | `ANT(2'')-Ia` = `aadB` |

**The `debla` rule** is the single highest-yield rule and reflects a real, systematic
nomenclature difference: CARD names beta-lactamases by allele only (`SHV-52`, `OXA-402`,
`CTX-M-58`), while ResFinder retains the `bla` gene-symbol prefix (`blaSHV-52`,
`blaOXA-402`, `blaCTX-M-58`) — both conventions are standard in the literature (Bush &
Jacoby 2010).

**The `van_structural` rule** was derived by enumerating all 61 `van`-family entries in
`card_prevalence.txt` and all 91 in `notes.txt` side by side: CARD spells out auxiliary
*vanRSHAXYZ* cluster genes descriptively (`"vanH gene in vanA cluster"`), ResFinder uses
a compact `VanH-A` form. The regex-based rule reproduces this correspondence exactly for
every cluster gene present in both files (vanA/B/D/E/F/G/M/N/O/P clusters).

**The `alias_bridge` rule** mines `phenotypes.txt`'s free-text `Notes` column, which
sometimes carries a curator-supplied `"Alternative name X, Y"` annotation (107 rows).
This surfaces historical gene names that predate systematic aminoglycoside-modifying-
enzyme (AME) nomenclature — most importantly `strA` = `APH(3'')-Ib`, `strB` = `APH(6)-Id`,
`aadA` = `ANT(3'')-Ia`, `aadB` = `ANT(2'')-Ia` — a correspondence independently documented
in the foundational AME nomenclature literature (Shaw et al. 1993; Ramirez & Tolmasky
2010).

## 3. Two real errors caught during validation (not just typos — worth flagging explicitly)

**a) Prime-mark collapse.** An early version of the punctuation-normalization rule
stripped apostrophes entirely, which wrongly collapsed `APH(3')-Ia` and `APH(3'')-Ia`
into the same key. These are *different enzymes*: prime-mark count denotes hydroxyl
regiospecificity in AME nomenclature (Shaw et al. 1993) and is not cosmetic punctuation.
Fixed by preserving apostrophe count in the normalization key.

**b) Cross-class false positive.** CARD's `FAR-1` (ARO Categories: *antibiotic
inactivation; penam* — a beta-lactamase) string-collided with ResFinder's `far1` (an
unrelated *Fusidic acid resistance*-class gene) under case/punctuation-insensitive
matching. Caught by a **data-driven drug-class consistency check**: the pipeline learns
which CARD `ARO Categories` terms co-occur with which ResFinder `notes_class` values from
already-unambiguous high-confidence matches, then flags candidates that contradict that
learned crosswalk. `far1` was correctly excluded; the true match, `blaFAR-1`, was
promoted. The same check correctly demoted CARD's `SatA` (categories: *nucleoside
antibiotic*) to "no reliable match" rather than accepting its only alias-bridge candidate,
ResFinder's `vat(D)` (a streptogramin A acetyltransferase) — a plausible but unconfirmed
pairing that the class mismatch correctly flagged as untrustworthy.

Both cases, plus all 21 CARD entries where the raw candidate pool had more than one
member (mostly legitimate documented synonym pairs, e.g. **KPC-1/KPC-2** and
**OXA-23/ARI-1**-style historical renamings — see `ambiguous_cases_review.csv`), are
listed with full provenance for manual audit.

## 4. Validation performed

- Random sample of 6 matches drawn from **every** tier, checked by hand: 36/36 correct.
- All 140 chromosomal/regulatory-model CARD entries confirmed to produce zero matches.
- Unmatched `protein homolog model` entries (1,547) profiled: 491 attributable to
  species outside ResFinder's panels (PDC/ADC/OXY/LEN/MIR/ACT families), a further ~49
  are CARD's species-qualified efflux/regulatory genes (e.g. `MexC`, `Escherichia coli
  mdfA`) that fall outside ResFinder's acquired-gene, BLAST-based detection model
  entirely (efflux pumps are multi-component chromosomal systems, not single
  horizontally-transferred genes).
- One manual, literature-documented curation override applied (`APH(3'')-Ib` &rarr;
  `strA`, not `aph(6)-Ia`; both are only reachable via a shared third-party alias in
  ResFinder's own notes, and only one is the correct regiospecificity match) — logged
  explicitly in the script output.

## 5. Known limitations

- `card_prevalence.txt` is CARD's *NCBI-surveillance-detected* subset (ARO terms with
  &ge;1 perfect/perfect_strict RGI hit against NCBI plasmid/WGS/chromosome/genomic-island
  sequences), not the full CARD ontology. A ResFinder gene with no CARD counterpart here
  may still exist in full CARD but simply never registered a hit in the surveyed NCBI
  sequences.
- `alias_bridge` matches (Medium confidence) depend on the accuracy of ResFinder's own
  curator-supplied cross-references, which are occasionally internally inconsistent (see
  the `APH(3'')-Ib` case above) — flagged, not silently trusted.
- Remaining unmatched entries were not individually hand-verified beyond the categorical
  profiling in Section 4; a small number may reflect genuine naming conventions the
  automated rules don't yet cover, rather than true absence from ResFinder.

## 6. References

1. Alcock BP et al. **CARD 2023: expanded curation, support for machine learning, and
   resistome prediction at the Comprehensive Antibiotic Resistance Database.** *Nucleic
   Acids Res.* 2023;51(D1):D690-D699. doi:10.1093/nar/gkac920
2. Alcock BP et al. **CARD 2020: antibiotic resistome surveillance with the
   Comprehensive Antibiotic Resistance Database.** *Nucleic Acids Res.*
   2020;48(D1):D517-D525. doi:10.1093/nar/gkz935
3. Jia B et al. **CARD 2017: expansion and model-centric curation of the Comprehensive
   Antibiotic Resistance Database.** *Nucleic Acids Res.* 2017;45(D1):D566-D573.
4. McArthur AG et al. **The Comprehensive Antibiotic Resistance Database.** *Antimicrob
   Agents Chemother.* 2013;57(7):3348-3357.
5. Zankari E et al. **Identification of acquired antimicrobial resistance genes.**
   *J Antimicrob Chemother.* 2012;67(11):2640-2644. doi:10.1093/jac/dks261
6. Bortolaia V et al. **ResFinder 4.0 for predictions of phenotypes from genotypes.**
   *J Antimicrob Chemother.* 2020;75(12):3491-3500. doi:10.1093/jac/dkaa345
7. Zankari E et al. **PointFinder: a novel web tool for WGS-based detection of
   antimicrobial resistance associated with chromosomal point mutations in bacterial
   pathogens.** *J Antimicrob Chemother.* 2017;72(10):2764-2768. doi:10.1093/jac/dkx217
8. Shaw KJ, Rather PN, Hare RS, Miller GH. **Molecular genetics of aminoglycoside
   resistance genes and familial relationships of the aminoglycoside-modifying
   enzymes.** *Microbiol Rev.* 1993;57(1):138-163.
9. Ramirez MS, Tolmasky ME. **Aminoglycoside modifying enzymes.** *Drug Resist Updat.*
   2010;13(6):151-171.
10. Bush K, Jacoby GA. **An updated functional classification of beta-lactamases.**
    *Antimicrob Agents Chemother.* 2010;54(3):969-976.
11. Jacoby GA, Chow N, Waites KB. **Prevalence of plasmid-mediated quinolone
    resistance.** *Antimicrob Agents Chemother.* 2003;47(2):559-562. (qnr discovery /
    early nomenclature)
12. Liu YY et al. **Emergence of plasmid-mediated colistin resistance mechanism MCR-1 in
    animals and human beings in China: a microbiological and molecular biological
    study.** *Lancet Infect Dis.* 2016;16(2):161-168.
13. Feldgarden M et al. **AMRFinderPlus and the Reference Gene Catalog facilitate
    examination of the genomic links among antimicrobial resistance, stress response,
    and virulence.** *Sci Rep.* 2021;11:12728. (precedent for cross-database AMR gene
    reconciliation)

---
*Generated as part of an interactive CARD/ResFinder reconciliation session. See
`build_card_resfinder_mapping.py` for the complete, re-runnable pipeline.*

# Dowloading and updating NCBI BLAST databases.

---
## Introduction

NCBI databases are hosted online in the NCBI database.
URL: https://ftp.ncbi.nlm.nih.gov/blast/db/

A full description of its contents are found in: https://ftp.ncbi.nlm.nih.gov/blast/db/README

---
## List NCBI databases

On the command line:
1. The list of accessible databases can be accessed through:

```
module load blast/2.14.1+
update_blastdb.pl --showall
```
---
## Download NCBI‑nt (nucleotide database)

```
update_blastdb.pl --decompress nt
```
This downloads:
 - nt.xx.tar.gz files
 - decompresses them
 - builds the final BLAST DB folder with v5 format
 - includes taxonomic metadata automatically

Result:
A folder with nt.nsq, nt.nin, nt.ndb, nt.nhr, etc.

## Downloading alternative NCBI databases:
 - core_nt: 
```
update_blastdb.pl --decompress core_nt
```
 - refseq_rna:
```
update_blastdb.pl --decompress refseq_rna
```

## Summary of Major NCBI BLAST Databases

| Database                    | Database Type                  | Contents                                                                                              | Typical Use                                               | Pros                                                                            | Cons                                                                                |
| --------------------------- | ------------------------------ | ----------------------------------------------------------------------------------------------------- | --------------------------------------------------------- | ------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| **nt**                      | Nucleotide                     | Comprehensive nucleotide sequences from GenBank, RefSeq, EMBL, DDBJ, patents, environmental sequences | Taxonomic identification, genome annotation, metagenomics | Most comprehensive nucleotide resource; highest sensitivity for homology search | Very large (>100 GB); redundant; contains low-quality or poorly annotated sequences |
| **core_nt**                 | Nucleotide                     | Curated subset of `nt` containing higher-quality sequences                                            | Faster nucleotide homology searches                       | Smaller and faster than `nt`; fewer low-quality records                         | Slightly reduced sensitivity; some sequences missing                                |
| **nt_prok**                 | Nucleotide                     | Subset of `nt` restricted to bacteria and archaea                                                     | Microbial genome annotation                               | Faster microbial searches; relevant taxonomic space                             | Limited to prokaryotes                                                              |
| **nt_euk**                  | Nucleotide                     | Subset of `nt` restricted to eukaryotes                                                               | Eukaryotic genome or transcript searches                  | Faster than `nt`; eukaryote-focused                                             | Still large; redundancy remains                                                     |
| **nt_viruses**              | Nucleotide                     | Viral nucleotide sequences from `nt`                                                                  | Viral detection, pathogen identification                  | Focused viral database; smaller than nt                                         | Misses unannotated viral sequences                                                  |
| **env_nt**                  | Nucleotide                     | Environmental and metagenomic nucleotide sequences                                                    | Environmental microbiology, novel diversity discovery     | Contains uncultured/environmental diversity                                     | Often poorly annotated; noisy taxonomy                                              |
| **nr**                      | Protein                        | Non-redundant protein sequences from GenBank CDS translations, RefSeq proteins, PDB, SwissProt, PIR   | Protein homology, functional annotation                   | Largest protein database; high sensitivity for distant homologs                 | Extremely large; annotation inconsistencies; redundancy by sequence clusters        |
| **env_nr**                  | Protein                        | Environmental metagenomic protein sequences                                                           | Discovery of novel proteins                               | Contains uncultured protein diversity                                           | Sparse annotation; noisy hits                                                       |
| **refseq_protein**          | Protein                        | Curated protein sequences from RefSeq genomes                                                         | High-quality protein annotation                           | Consistent annotation; less redundancy than nr                                  | Smaller database; may miss novel proteins                                           |
| **refseq_rna**              | Nucleotide                     | Curated RNA sequences from RefSeq                                                                     | Transcript annotation                                     | High-quality reference transcripts                                              | Limited diversity relative to nt                                                    |
| **refseq_select_prot**      | Protein                        | Representative curated protein per gene                                                               | Fast annotation pipelines                                 | Small, highly curated; consistent gene representation                           | Reduced coverage of sequence diversity                                              |
| **refseq_select_rna**       | Nucleotide                     | Representative transcript sequences                                                                   | Fast transcript annotation                                | Small, well curated                                                             | Limited isoform diversity                                                           |
| **ref_prok_rep_genomes**    | Genome (nucleotide + proteins) | Representative bacterial and archaeal genomes                                                         | Microbial genome comparison                               | One genome per species; curated                                                 | Limited strain diversity                                                            |
| **ref_euk_rep_genomes**     | Genome                         | Representative eukaryotic genomes                                                                     | Comparative genomics                                      | High-quality assemblies                                                         | Limited species coverage                                                            |
| **ref_viruses_rep_genomes** | Genome                         | Curated viral reference genomes                                                                       | Virus detection pipelines                                 | Reliable references                                                             | Misses many environmental viruses                                                   |
| **swissprot**               | Protein                        | Manually curated protein sequences from UniProtKB/Swiss-Prot                                          | Functional annotation                                     | Highest annotation quality; curated functions                                   | Small relative to nr                                                                |
| **pdbaa**                   | Protein                        | Protein sequences corresponding to structures in the Protein Data Bank                                | Structure-based analysis                                  | Links sequence to 3D structure                                                  | Limited to proteins with solved structures                                          |
| **pdbnt**                   | Nucleotide                     | Nucleotide sequences associated with PDB entries                                                      | Structural genomics                                       | Direct link to structural data                                                  | Very small database                                                                 |
| **tsa_nr**                  | Protein                        | Translated proteins from transcriptome shotgun assemblies                                             | Non-model organism transcriptomics                        | Expands coverage of species                                                     | Annotation often unreliable                                                         |
| **tsa_nt**                  | Nucleotide                     | Transcriptome shotgun assemblies                                                                      | Transcript discovery                                      | Large coverage of transcriptomes                                                | Fragmented sequences                                                                |


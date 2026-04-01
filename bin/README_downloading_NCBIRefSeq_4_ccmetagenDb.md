# **How to Download and Format NCBI's RefSeq for CCMETAGEN DB**

This has been executed in the script **`download_all_refseq.sh`**.  
The script automates the creation of a **RefSeq‑based database** suitable for **CCMetagen**, **metagenomic classification**, and **taxonomic annotation**. It performs controlled downloading, retry handling, taxonomic mapping, and header normalization to produce a clean, taxid‑labeled FASTA reference.

---

## **1. Download the RefSeq assembly summary**

The script retrieves (or reuses) the master RefSeq `assembly_summary_refseq.txt`, which describes all assemblies in the NCBI RefSeq database.  
This file is the **source of all genome metadata**, including which genomes are labeled *reference* and where to find their FTP paths.

---

## **2. Extract and clean the list of genomes to download**

Using column‑based filtering, the script isolates only:

*   **Reference genomes** (*and optionally representative genomes*)
*   Their **FTP paths** (column 20)

The URLs are cleaned (quotes and trailing slashes removed), producing a **canonical list** of RefSeq genome directories to fetch.

---

## **3. Polite, fault‑tolerant multi‑round genome downloader**

Genomes are downloaded in a **controlled**, **polite**, **retry‑aware** process:

*   **Parallelized** (default 8 jobs) with short delays to avoid hammering NCBI.
*   Up to **10 rounds** of retries, with **exponential backoff**.
*   Each genome:
    *   Is **skipped** if already downloaded.
    *   Is fetched to a temporary file, then atomically moved into place.
    *   Produces clear `STATUS:` logs (`OK`, `GET`, `DONE`, `FAIL`).
*   MD5 checksum entries from each genome directory are captured for reference.

This stage ensures high reliability even under network congestion or NCBI throttling.

---

## **4. Download NCBI accession→taxid mapping files**

The script retrieves the massive `nucl_gb` and `nucl_wgs` accession‑to‑taxid maps.  
These contain authoritative links between NCBI accessions and taxonomy IDs—critical for labeling genomes correctly.

---

## **5. Map genome FASTAs to their accession numbers**

For each downloaded genome FASTA:

*   The script extracts the **accession** from its header.
*   Creates a table linking **accession → file path**.

This allows joining taxonomy metadata to the downloaded files.

---

## **6. Build a unified taxid→file mapping**

Using a join between:

*   **accession→taxid** map (from NCBI metadata)
*   **accession→file** map (from downloaded FASTAs)

…the script produces a clean mapping:

    taxid    path/to/genome.fna.gz

This file drives the header‑restamping step.

---

## **7. Restamp FASTA headers into CCMetagen format**

Each genome FASTA is transformed so that every sequence header becomes:

    >taxid|original_accession_and_description

This is the **required CCMetagen format**, enabling fast and accurate taxonomic classification.

The script:

*   Streams each FASTA through `sed` to rewrite headers.
*   Writes new files ending in `_taxID2.fna`.
*   Compresses the results and removes the original FASTAs.

---

## **8. Final deliverables**

At the end, the directory contains:

*   **Taxid‑labeled FASTA files** (`*_taxID2.fna.gz`)
*   **MD5 metadata** per download round
*   **Accession→taxid**, **accession→file**, and **taxid→file** maps
*   **Download logs**, **parallel logs**, and **failure lists** (if any)

These files collectively form a **clean, validated, taxonomically annotated RefSeq reference database**.

---

# **Why these steps matter**

*   Ensures **reliable downloading** from NCBI with retries and throttling‑safe behavior.
*   Guarantees **traceability** through MD5 checks and comprehensive logs.
*   Produces **consistent taxonomic labeling**, essential for CCMetagen and other classifiers.
*   Creates a **restartable workflow**, safe against partial downloads or interruptions.
*   Outputs a fully normalized, taxonomically annotated genome database.

---

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/ccmetagen_db/main.nf

    Prepare a CCMetagen / KMA database from one of four sources:

      PRIORITY 1 — ZIP archive  (local .zip file)
          → unzips into db_out/

      PRIORITY 2 — Pre-built KMA DB directory  (existing indexed KMA dir)
          → copies into db_out/ without rebuilding

      PRIORITY 3 — Pre-cleaned FASTA file(s)  ← NEW: your local use case
          → runs `kma index` directly on .fasta.gz / .fna.gz / .fa.gz
          → use this when you already have a cleaned, taxonomy-annotated FASTA
            (e.g. core_nt.taxid_desc.ccmetagen.fasta.gz) and just need indexing

      PRIORITY 4a — Keyword + downloadccmetagen_db
          → wget pre-built archive from Unimelb mediaflux + unzip

      PRIORITY 4b — Keyword + buildccmetagen_db
          → download raw FASTA from NCBI/RefSeq then kma index

    ── Input ───────────────────────────────────────────────────────────────────
    tuple val(db_id), path(db_in)
        db_in = one of:
          • path/to/archive.zip            → PRIORITY 1
          • path/to/existing_kma_dir/      → PRIORITY 2
          • path/to/cleaned.fasta.gz       → PRIORITY 3  ← your case
          • file("NO_FILE")                → PRIORITY 4a / 4b (keyword mode)

    val  db_action    // null | "downloadccmetagen_db" | "buildccmetagen_db"
    val  db_keyword   // only for PRIORITY 4: "RefSeq" | "NCBI_nt" | "NCBI_nt_no_env"

    ── Output ──────────────────────────────────────────────────────────────────
    tuple val(db_id), path("db_out")   , emit: db
    path  "versions.yml"               , emit: versions

    ── Quick usage examples ─────────────────────────────────────────────────────
    // PRIORITY 3 — index a pre-cleaned FASTA (your current situation):
    CCMETAGEN_DB(
        tuple("core_nt", file("ccmetagen_db/core_nt/core_nt.taxid_desc.ccmetagen.fasta.gz")),
        null,
        null
    )

    // PRIORITY 2 — reuse an already-indexed KMA directory:
    CCMETAGEN_DB(
        tuple("core_nt", file("path/to/kma_db_dir/")),
        null,
        null
    )

    // PRIORITY 4a — download the Unimelb pre-built archive:
    CCMETAGEN_DB(
        tuple("NCBI_nt", file("NO_FILE")),
        "downloadccmetagen_db",
        "NCBI_nt"
    )
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CCMETAGEN_DB {

    tag "${db_id}"
    label 'process_high_memory'

    // KMA container — used for PRIORITY 2, 3, 4b.
    // PRIORITY 1 and 4a also need wget+unzip; those tools are available in the
    // base image below. Uncomment one profile once you confirm which tool stack
    // your HPC uses.
    conda "bioconda::kma=1.4.14 conda-forge::wget=1.21"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.14--h9f5acd7_0' :
        'quay.io/biocontainers/kma:1.4.14--h9f5acd7_0' }"

    // Persist the final DB directory — critical for an expensive build step.
    publishDir(
        path:    { "${params.outdir}/ccmetagen_db/${db_id}" },
        mode:    params.publish_dir_mode ?: 'copy',
        pattern: "db_out/**"
    )

    input:
    tuple val(db_id), path(db_in)   // ZIP, KMA-dir, FASTA, or NO_FILE placeholder
    val   db_action                 // null | "downloadccmetagen_db" | "buildccmetagen_db"
    val   db_keyword                // "RefSeq" | "NCBI_nt" | "NCBI_nt_no_env" (PRIORITY 4 only)

    output:
    tuple val(db_id), path("db_out"), emit: db
    path  "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def taxon      = params.ccmetagen_taxon ?: 'bacteria'
    // kma index is single-threaded so -t is NOT forwarded to the index command.
    // Retained here for potential use in build mode's download/preprocessing steps.
    def kma_threads = task.cpus

    // ── Remote URLs ──────────────────────────────────────────────────────────
    def db_urls = [
        "RefSeq"        : "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=Lqaic1pBmpDdqX8ofv1C1128247855&browser=true&filename=RefSeq_bf.zip",
        "NCBI_nt"       : "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=i8yedNiYfdjrBfGJ8Y5z1128247857&browser=true&filename=ncbi_nt_kma.zip",
        "NCBI_nt_no_env": "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=ko6MbZXl7FWjAS3jsItV1128247851&browser=true&filename=ncbi_nt_no_env_11jun2019.zip"
    ]

    // ── NCBI FASTA sources for build-from-scratch ────────────────────────────
    def fasta_sources = [
        "NCBI_nt": "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz"
    ]

    // ── Input classification ─────────────────────────────────────────────────
    def db_name  = db_in.name
    def is_zip   = db_name.endsWith(".zip")
    def is_dir   = db_in.isDirectory()

    // Detect pre-cleaned FASTA files so they route to PRIORITY 3
    // instead of falling through to the else-error block.
    def fasta_extensions = ['.fasta.gz', '.fna.gz', '.fa.gz', '.fasta', '.fna', '.fa']
    def is_fasta  = !is_zip && !is_dir && fasta_extensions.any { db_name.endsWith(it) }

    // is_keyword only when none of the above matched
    def is_keyword = (!is_zip && !is_dir && !is_fasta && db_keyword)

    // ── PRIORITY 1: ZIP archive ───────────────────────────────────────────────
    if ( is_zip ) {
        """
        set -euo pipefail
        echo "[CCMETAGEN_DB] PRIORITY 1 — Unzipping archive: ${db_in}"
        mkdir -p db_out
        unzip -q "${db_in}" -d db_out/
        echo "[CCMETAGEN_DB] Unzip complete."

        cat <<END_VERSIONS > versions.yml
        "${task.process}":
            unzip: \$(unzip -v 2>&1 | head -1 | awk '{print \$2}')
        END_VERSIONS
        """

    // ── PRIORITY 2: existing pre-built KMA DB directory ──────────────────────
    } else if ( is_dir ) {
        """
        set -euo pipefail
        echo "[CCMETAGEN_DB] PRIORITY 2 — Using pre-built KMA DB directory: ${db_in}"
        mkdir -p db_out
        # cp -rL dereferences symlinks so the staged copy is fully self-contained
        cp -rL "${db_in}/." db_out/
        echo "[CCMETAGEN_DB] Copy complete. Contents:"
        ls -lh db_out/

        cat <<END_VERSIONS > versions.yml
        "${task.process}":
            kma: \$(kma -v 2>&1 | head -1 | sed 's/KMA-//' || echo "N/A")
        END_VERSIONS
        """

    // ── PRIORITY 3: pre-cleaned FASTA → kma index ────────────────────────────
    // Use this when you already have a taxonomy-annotated, deduplicated FASTA
    // and just need it indexed for CCMetagen classification.
    // The db_id is used as the index prefix so KMA output files are named
    // consistently (e.g. db_out/core_nt.name, db_out/core_nt.seq, …).
    } else if ( is_fasta ) {
        """
        set -euo pipefail
        echo "[CCMETAGEN_DB] PRIORITY 3 — Indexing pre-cleaned FASTA: ${db_in}"
        echo "[CCMETAGEN_DB] Output prefix  : db_out/${db_id}"
        echo "[CCMETAGEN_DB] KMA version    : \$(kma -v 2>&1 | head -1)"

        mkdir -p db_out

        kma index \\
            -i "${db_in}" \\
            -o db_out/${db_id} \\
            ${args}

        echo "[CCMETAGEN_DB] Indexing complete. Index files:"
        ls -lh db_out/

        cat <<END_VERSIONS > versions.yml
        "${task.process}":
            kma: \$(kma -v 2>&1 | head -1 | sed 's/KMA-//')
        END_VERSIONS
        """

    // ── PRIORITY 4a: download pre-built archive by keyword ───────────────────
    } else if ( is_keyword && db_action == "downloadccmetagen_db" ) {
        def url = db_urls[db_keyword]
        if ( !url ) {
            error "[CCMETAGEN_DB] Unknown keyword for download: '${db_keyword}'. Valid: ${db_urls.keySet()}"
        }

        """
        set -euo pipefail
        echo "[CCMETAGEN_DB] PRIORITY 4a — Downloading pre-built CCMetagen DB: ${db_keyword}"
        wget -q --show-progress -O db.zip "${url}"
        mkdir -p db_out
        unzip -q db.zip -d db_out/
        rm -f db.zip
        echo "[CCMETAGEN_DB] Download and unzip complete."

        cat <<END_VERSIONS > versions.yml
        "${task.process}":
            wget: \$(wget --version 2>&1 | head -1 | sed 's/GNU Wget //')
            unzip: \$(unzip -v 2>&1 | head -1 | awk '{print \$2}')
        END_VERSIONS
        """

    // ── PRIORITY 4b: build from scratch — RefSeq ─────────────────────────────
    } else if ( is_keyword && db_action == "buildccmetagen_db" && db_keyword == "RefSeq" ) {

        """
        set -euo pipefail
        echo "[CCMETAGEN_DB] PRIORITY 4b — Building CCMetagen RefSeq DB for taxon: ${taxon}"
        mkdir -p db_out

        download_all_refseq.sh "${taxon}"

        echo "[CCMETAGEN_DB] Concatenating downloaded FASTA files..."
        find . -name "*.fna.gz" -exec gunzip -c {} \\; > refseq_all.fna

        echo "[CCMETAGEN_DB] Running kma index..."
        kma index \\
            -i refseq_all.fna \\
            -o db_out/${db_keyword} \\
            ${args}

        rm -f refseq_all.fna
        echo "[CCMETAGEN_DB] RefSeq database build complete."

        cat <<END_VERSIONS > versions.yml
        "${task.process}":
            kma: \$(kma -v 2>&1 | head -1 | sed 's/KMA-//')
        END_VERSIONS
        """

    // ── PRIORITY 4b: build from scratch — NCBI_nt or custom FASTA URL ────────
    } else if ( is_keyword && db_action == "buildccmetagen_db" ) {
        def fasta = fasta_sources[db_keyword]
        if ( !fasta ) {
            error "[CCMETAGEN_DB] No FASTA source for keyword '${db_keyword}'. Available: ${fasta_sources.keySet()}"
        }

        """
        set -euo pipefail
        echo "[CCMETAGEN_DB] PRIORITY 4b — Building KMA DB for: ${db_keyword}"
        mkdir -p db_out

        echo "[CCMETAGEN_DB] Downloading FASTA: ${fasta}"
        wget -q --show-progress "${fasta}" -O source.fna.gz

        echo "[CCMETAGEN_DB] Running kma index..."
        kma index \\
            -i source.fna.gz \\
            -o db_out/${db_keyword} \\
            ${args}

        rm -f source.fna.gz
        echo "[CCMETAGEN_DB] Database build complete."

        cat <<END_VERSIONS > versions.yml
        "${task.process}":
            kma: \$(kma -v 2>&1 | head -1 | sed 's/KMA-//')
            wget: \$(wget --version 2>&1 | head -1 | sed 's/GNU Wget //')
        END_VERSIONS
        """

    } else {
        error """
        [CCMETAGEN_DB] Unsupported DB configuration:
          db_id      = ${db_id}
          db_in      = ${db_in}
          db_action  = ${db_action}
          db_keyword = ${db_keyword}

        Valid combinations:
          (1) db_in = path/to/archive.zip                    → unzip
          (2) db_in = path/to/existing_kma_dir/              → use as-is
          (3) db_in = path/to/cleaned.fasta.gz               → kma index   ← your case
          (4a) db_keyword + db_action=downloadccmetagen_db   → wget prebuilt
          (4b) db_keyword + db_action=buildccmetagen_db      → download FASTA + kma index
        """
    }

    stub:
    """
    mkdir -p db_out
    touch db_out/${db_id}.name db_out/${db_id}.seq db_out/${db_id}.len db_out/${db_id}.comp.b
    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        kma: 1.4.14
    END_VERSIONS
    """
}

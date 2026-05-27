/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/ccmetagen_db.nf
    Prepare a CCMetagen / KMA database from one of three sources:

      PRIORITY 1 — ZIP archive  (local file ending in .zip)
          → unzips into db_out/
      PRIORITY 2 — Pre-built directory  (path to existing KMA DB dir)
          → symlinks/copies into db_out/ without rebuilding
      PRIORITY 3 — Keyword + action
          a. downloadccmetagen_db  → wget pre-built archive from Unimelb mediaflux
          b. buildccmetagen_db     → build from FASTA (RefSeq download or NCBI_nt)

    Input:
        tuple val(db_id), path(db_in)   // zip/dir path, OR a dummy file when keyword
        val   db_action                 // null | "downloadccmetagen_db" | "buildccmetagen_db"
    Output:
        tuple val(db_id), path("db_out")    , emit: db
        path  "versions.yml"                , emit: versions

    Note: for keyword-mode (db_in is a dummy file), pass:
        db_in     = file("NO_FILE")  // placeholder
        db_action = "downloadccmetagen_db" | "buildccmetagen_db"
        And set  params.ccmetagen_db  to the keyword string (RefSeq | NCBI_nt | NCBI_nt_no_env)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CCMETAGEN_DB {

    tag "${db_id}"

    label 'process_high_memory'

    // conda "bioconda::ccmetagen=1.5.0 conda-forge::wget=1.21"
    // container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    //     'oras://community.wave.seqera.io/library/kma:1.4.14' :
    //     'biocontainers/kma:1.4.14' }"

    input:
    tuple val(db_id), path(db_in)    // ZIP file, existing DB directory, or NO_FILE placeholder
    val   db_action                  // null | "downloadccmetagen_db" | "buildccmetagen_db"
    val   db_keyword                 // keyword string when db_in is NO_FILE: RefSeq|NCBI_nt|…

    output:
    tuple val(db_id), path("db_out") , emit: db
    path  "versions.yml"             , emit: versions

    script:
    def args       = task.ext.args   ?: ''
    def taxon      = params.ccmetagen_taxon ?: 'bacteria'
    def kma_threads = task.cpus

    // Remote pre-built CCMetagen DB URLs (Unimelb mediaflux)
    def db_urls = [
        "RefSeq"        : "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=Lqaic1pBmpDdqX8ofv1C1128247855&browser=true&filename=RefSeq_bf.zip",
        "NCBI_nt"       : "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=i8yedNiYfdjrBfGJ8Y5z1128247857&browser=true&filename=ncbi_nt_kma.zip",
        "NCBI_nt_no_env": "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=ko6MbZXl7FWjAS3jsItV1128247851&browser=true&filename=ncbi_nt_no_env_11jun2019.zip"
    ]

    // FASTA sources for build-from-scratch
    def fasta_sources = [
        "NCBI_nt" : "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz"
    ]

    def db_name    = db_in.name     // "NO_FILE" when keyword mode
    def is_zip     = db_name.endsWith(".zip")
    def is_dir     = db_in.isDirectory()
    def is_keyword = (!is_zip && !is_dir && db_keyword)

    if ( is_zip ) {
        // ── PRIORITY 1: ZIP ────────────────────────────────────────────────
        """
        echo "[CCMETAGEN_DB] Unzipping provided archive: ${db_in}"
        mkdir -p db_out
        unzip -q "${db_in}" -d db_out/
        echo "Done."

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            unzip: \$(unzip -v 2>&1 | head -1 | awk '{print \$2}')
        END_VERSIONS
        """

    } else if ( is_dir ) {
        // ── PRIORITY 2: existing pre-built directory ───────────────────────
        """
        echo "[CCMETAGEN_DB] Using pre-built DB directory: ${db_in}"
        # Create db_out as a symlink-copy so Nextflow stages it correctly
        mkdir -p db_out
        cp -rL "${db_in}/." db_out/
        echo "Done."

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kma: \$(kma -v 2>&1 | head -1 | sed 's/KMA-//' || echo "N/A")
        END_VERSIONS
        """

    } else if ( is_keyword && db_action == "downloadccmetagen_db" ) {
        // ── PRIORITY 3a: download pre-built archive by keyword ─────────────
        def url = db_urls[db_keyword]
        if ( !url ) error "[CCMETAGEN_DB] Unknown keyword for download: '${db_keyword}'. Valid: ${db_urls.keySet()}"

        """
        echo "[CCMETAGEN_DB] Downloading pre-built CCMetagen DB for '${db_keyword}'"
        wget -q --show-progress -O db.zip "${url}"
        mkdir -p db_out
        unzip -q db.zip -d db_out/
        rm -f db.zip
        echo "Done."

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            wget: \$(wget --version 2>&1 | head -1 | sed 's/GNU Wget //')
            unzip: \$(unzip -v 2>&1 | head -1 | awk '{print \$2}')
        END_VERSIONS
        """

    } else if ( is_keyword && db_action == "buildccmetagen_db" ) {
        // ── PRIORITY 3b: build DB from scratch via KMA index ──────────────

        if ( db_keyword == "RefSeq" ) {
            // Download RefSeq genomes per taxon then KMA-index
            """
            echo "[CCMETAGEN_DB] Building CCMetagen RefSeq DB for taxon: ${taxon}"
            mkdir -p db_out

            # Download RefSeq genomes (uses bin/download_all_refseq.sh)
            download_all_refseq.sh "${taxon}"

            echo "Concatenating downloaded FASTA files..."
            find . -name "*.fna.gz" -exec gunzip -c {} \\; > refseq_all.fna

            echo "Running KMA index..."
            kma index \\
                -i refseq_all.fna \\
                -o db_out/${db_keyword} \\
                -t_db db_out \\
                -t ${kma_threads} \\
                ${args}

            rm -f refseq_all.fna
            echo "Database build complete."

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                kma: \$(kma -v 2>&1 | head -1 | sed 's/KMA-//')
            END_VERSIONS
            """

        } else {
            // NCBI_nt or similar — KMA-index from FASTA URL
            def fasta = fasta_sources[db_keyword]
            if ( !fasta ) error "[CCMETAGEN_DB] No FASTA source for keyword '${db_keyword}'. Available: ${fasta_sources.keySet()}"

            """
            echo "[CCMETAGEN_DB] Building CCMetagen/KMA DB for '${db_keyword}'"
            mkdir -p db_out

            echo "Downloading FASTA: ${fasta}"
            wget -q --show-progress "${fasta}" -O source.fna.gz

            echo "Running KMA index..."
            kma index \\
                -i source.fna.gz \\
                -o db_out/${db_keyword} \\
                -t_db db_out \\
                -t ${kma_threads} \\
                ${args}

            rm -f source.fna.gz
            echo "Database build complete."

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                kma: \$(kma -v 2>&1 | head -1 | sed 's/KMA-//')
                wget: \$(wget --version 2>&1 | head -1 | sed 's/GNU Wget //')
            END_VERSIONS
            """
        }

    } else {
        error """
        [CCMETAGEN_DB] Unsupported DB configuration:
          db_id      = ${db_id}
          db_in      = ${db_in}
          db_action  = ${db_action}
          db_keyword = ${db_keyword}

        Valid combinations:
          (a) db_in = path/to/file.zip          → unzip
          (b) db_in = path/to/existing_db_dir   → use as-is
          (c) db_keyword + db_action = "downloadccmetagen_db"   → wget prebuilt
          (d) db_keyword + db_action = "buildccmetagen_db"      → build from FASTA
        """
    }

    stub:
    """
    mkdir -p db_out
    touch db_out/stub.name db_out/stub.seq db_out/stub.len
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: 1.4.14
    END_VERSIONS
    """
}

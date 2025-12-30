process CCMETAGEN_DB {

    tag "${db_id}"

    input:
    tuple val(db_id), path(db_in) // directory, zip, or keyword (e.g. "RefSeq")
    val db_action          // null, "buildccmetagen_db", "downloadccmetagen_db"

    output:
    tuple val(db_id), path("db_out")

    script:

    // Detect input types according to required priority
    def isDir   = db_in instanceof Path && db_in.isDirectory()
    def isZip   = db_in instanceof Path && db_in.name.endsWith(".zip")
    def isKey   = db_in instanceof String && !db_in.contains("/") && !db_in.contains(".")

    // Remote prebuilt CCMetagen DB archives
    def db_urls = [
        "RefSeq"            : "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=Lqaic1pBmpDdqX8ofv1C1128247855&browser=true&filename=RefSeq_bf.zip",
        "NCBI_nt"           : "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=i8yedNiYfdjrBfGJ8Y5z1128247857&browser=true&filename=ncbi_nt_kma.zip",
        "NCBI_nt_no_env"    : "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=ko6MbZXl7FWjAS3jsItV1128247851&browser=true&filename=ncbi_nt_no_env_11jun2019.zip"
    ]

    // FASTA sources for building DBs from scratch
    def fasta_sources = [
        "NCBI_nt" : "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz"
    ]

    // -----------------------------------------------------------------
    // PRIORITY 1: ZIP file wins over everything else
    // -----------------------------------------------------------------
    if (isZip) {
        """
        echo "DB: Unzipping provided archive"
        mkdir db_out
        unzip -q ${db_in} -d db_out
        """
    }

    // -----------------------------------------------------------------
    // PRIORITY 2: Directory
    // -----------------------------------------------------------------
    else if (isDir) {
        """
        echo "DB: Using provided directory"
        ln -s ${db_in} db_out
        """
    }

    // -----------------------------------------------------------------
    // PRIORITY 3: Key + action
    // -----------------------------------------------------------------
    else if (isKey) {

        // Case A: download
        if (db_action == "downloadccmetagen_db") {
            def url = db_urls[db_in]
            if (!url) error "Unknown CCMetagen DB key for download: ${db_in}"

            """
            echo "DB: Downloading prebuilt CCMetagen DB for key '${db_in}'"
            wget -q -O db.zip "${url}"
            mkdir db_out
            unzip -q db.zip -d db_out
            """
        }

        // Case B: build (KMA indexing)
        else if (db_action == "buildccmetagen_db") {

            if (db_in == "RefSeq") {
                """
                echo "DB: Building CCMetagen RefSeq DB for taxon: ${params.ccmetagen_taxon}"
                mkdir db_out

                ../../bin/download_refseq.sh "${params.ccmetagen_taxon}"

                echo "Concatenating FASTA..."
                find refseq_*_fasta -name "*.fna.gz" -exec gunzip -c {} \; > refseq_all.fna

                kma index \
                    -i refseq_all.fna \
                    -o db_out/${db_in} \
                    -t_db db_out

                echo "Database build complete."
                """
            }

            else {
                def fasta = fasta_sources[db_in]
                if (!fasta) error "No FASTA source available for DB key '${db_in}'"

                """
                echo "DB: Building CCMetagen/KMA DB for key '${db_in}'"
                mkdir db_out

                kma index \
                    -i ${fasta} \
                    -o db_out/${db_in} \
                    -t_db db_out

                echo "Database build complete."
                """
            }
        }

        else {
            error "DB key provided but action is missing or invalid: db_in=${db_in}, db_action=${db_action}"
        }
    }

    // -----------------------------------------------------------------
    // NOT ZIP, NOT DIR, NOT KEY → error
    // -----------------------------------------------------------------
    else {
        error """
        Unsupported CCMetagen DB configuration:
        db_id     = ${db_id}
        db_in     = ${db_in}
        db_action = ${db_action}
        """
    }
}


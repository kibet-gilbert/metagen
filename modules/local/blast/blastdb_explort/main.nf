/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/blast_export.nf
    MODULE 1 of 2 in the CCMetagen DB-build chain.

    Exports raw sequences from an NCBI BLAST database (nt, core_nt, etc.)
    with CCMetagen-compatible headers ( >taxid|title ) using export_blastdb.sh.

    Output is a gzip-compressed FASTA ready for FASTA_CLEAN (module 2).

    Input:
        tuple val(meta), path(blast_db_dir)  // directory holding BLAST DB files
                                             // meta.id = short label e.g. "core_nt"
    Output:
        tuple val(meta), path("*.raw.fasta.gz")  , emit: raw_fasta
        path  "versions.yml"                     , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BLAST_EXPORT {

    tag "${meta.id}"

    label 'process_high'

    // Uses NCBI BLAST+ (blastdbcmd) and pigz — load via module or container.
    conda "bioconda::blast=2.15.0 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1' :
        'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta), path(blast_db_dir)   // BLAST DB directory (e.g. /data/ncbi/nt)

    output:
    tuple val(meta), path("*.raw.fasta.gz") , emit: raw_fasta
    path  "versions.yml"                    , emit: versions

    script:
    def args       = task.ext.args        ?: ''
    def prefix     = task.ext.prefix      ?: "${meta.id}"
    def compress   = task.ext.compress    ?: '1'
    def header_fmt = task.ext.header_fmt  ?: '>%T|%t\\n%s'

    // The script expects a path prefix (without extension), so we find
    // the actual DB prefix from the directory contents.
    """
    # ── locate the blast DB prefix inside the staged directory ─────────────
    DB_PREFIX=\$(find "${blast_db_dir}" \\
        -maxdepth 1 \\
        -name "*.nsi" -o -name "*.nhr" -o -name "*.nal" 2>/dev/null \\
        | head -1 \\
        | sed 's/\\.[^.]*\$//')

    if [[ -z "\$DB_PREFIX" ]]; then
        # Try passing the directory itself as the prefix (blastdbcmd handles it)
        DB_PREFIX="${blast_db_dir}/${prefix}"
    fi

    # ── run the export script ──────────────────────────────────────────────
    export HEADER_FMT='${header_fmt}'
    export THREADS="${task.cpus}"
    export COMPRESS="${compress}"

    export_blastdb.sh "\$DB_PREFIX" export_out/ ${args}

    # ── stage output back to work dir with consistent name ─────────────────
    mv export_out/*.raw.fasta.gz ./ 2>/dev/null \\
        || mv export_out/*.raw.fasta ./ 2>/dev/null \\
        || { echo "ERROR: no FASTA output found in export_out/"; exit 1; }

    # Ensure output is always gzipped
    for f in *.raw.fasta; do
        [[ -f "\$f" ]] && pigz -p "${task.cpus}" "\$f"
    done

    # ── versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastdbcmd: \$(blastdbcmd -version 2>&1 | head -1 | sed 's/blastdbcmd: //')
        pigz: \$(pigz --version 2>&1 | head -1 || echo "N/A")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "stub" | gzip > ${prefix}.raw.fasta.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastdbcmd: 2.15.0
        pigz: 2.8
    END_VERSIONS
    """
}

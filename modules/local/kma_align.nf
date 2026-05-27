/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/kma_align.nf
    Align metagenomic reads to a KMA database (CCMetagen-compatible).

    KMA outputs two key files per sample that CCMetagen.py consumes:
      *.res      — per-taxon read-count / coverage / depth table
      *.mapstat  — per-read mapping statistics

    Input:
        tuple val(meta), path(reads)     // [R1, R2] paired, or [reads] single
        tuple val(db_id), path(db_dir)   // from CCMETAGEN_DB or FASTA_CLEAN
    Output:
        tuple val(meta), path("*.res")                 , emit: res
        tuple val(meta), path("*.mapstat")             , emit: mapstat
        tuple val(meta), path("*.res"), path("*.mapstat"), emit: res_mapstat
        path  "versions.yml"                           , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KMA_ALIGN {

    tag "${meta.id}"

    label 'process_high'

    // conda "bioconda::kma=1.4.14"
    // container "${ workflow.containerEngine == 'singularity' ?
    //     'oras://community.wave.seqera.io/library/kma:1.4.14' :
    //     'biocontainers/kma:1.4.14' }"

    input:
    tuple val(meta), path(reads)
    tuple val(db_id), path(db_dir)

    output:
    tuple val(meta), path("${meta.id}.res")                          , emit: res
    tuple val(meta), path("${meta.id}.mapstat")                      , emit: mapstat
    tuple val(meta), path("${meta.id}.res"), path("${meta.id}.mapstat"), emit: res_mapstat
    path  "versions.yml"                                             , emit: versions

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"

    // Locate the KMA DB prefix inside the staged directory
    // The DB directory will contain files like: <name>.name, <name>.seq, <name>.len
    def paired = meta.single_end ? false : true

    if ( paired ) {
        """
        # ── Locate KMA DB prefix ──────────────────────────────────────────
        DB_PREFIX=\$(find "${db_dir}" \\
            -maxdepth 2 -name "*.name" 2>/dev/null \\
            | head -1 | sed 's/\\.name\$//')

        if [[ -z "\$DB_PREFIX" ]]; then
            echo "ERROR: Could not locate KMA DB prefix in ${db_dir}"
            ls -lh "${db_dir}/"
            exit 1
        fi

        echo "[KMA] DB prefix : \$DB_PREFIX"
        echo "[KMA] Input     : ${reads[0]} ${reads[1]}"

        kma \\
            -ipe ${reads[0]} ${reads[1]} \\
            -o ${prefix} \\
            -t_db "\$DB_PREFIX" \\
            -t ${task.cpus} \\
            -mem_mode \\
            -and \\
            -apm f \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kma: \$(kma -v 2>&1 | head -1 | sed 's/KMA-//')
        END_VERSIONS
        """
    } else {
        """
        # ── Locate KMA DB prefix ──────────────────────────────────────────
        DB_PREFIX=\$(find "${db_dir}" \\
            -maxdepth 2 -name "*.name" 2>/dev/null \\
            | head -1 | sed 's/\\.name\$//')

        if [[ -z "\$DB_PREFIX" ]]; then
            echo "ERROR: Could not locate KMA DB prefix in ${db_dir}"
            ls -lh "${db_dir}/"
            exit 1
        fi

        echo "[KMA] DB prefix : \$DB_PREFIX"
        echo "[KMA] Input     : ${reads}"

        kma \\
            -i ${reads} \\
            -o ${prefix} \\
            -t_db "\$DB_PREFIX" \\
            -t ${task.cpus} \\
            -mem_mode \\
            -and \\
            -apm f \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kma: \$(kma -v 2>&1 | head -1 | sed 's/KMA-//')
        END_VERSIONS
        """
    }

    stub:
    """
    touch ${meta.id}.res ${meta.id}.mapstat
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: 1.4.14
    END_VERSIONS
    """
}

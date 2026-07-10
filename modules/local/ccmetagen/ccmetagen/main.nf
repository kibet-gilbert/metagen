/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/ccmetagen.nf
    Run CCMetagen.py on per-sample KMA alignment results.

    CCMetagen.py takes the .res file (+ optional .mapstat) produced by KMA and
    produces a per-sample taxonomic abundance table in CSV format.

    Input:
        tuple val(meta), path(res), path(mapstat)   // from KMA_ALIGN
        tuple val(db_id), path(db_dir)              // from CCMETAGEN_DB / FASTA_CLEAN
    Output:
        tuple val(meta), path("${meta.id}_ccmetagen/")       , emit: results_dir
        tuple val(meta), path("${meta.id}_ccmetagen/*.csv")  , emit: csv
        path  "versions.yml"                                 , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CCMETAGEN {

    tag "${meta.id}"

    label 'process_medium'

    // conda "bioconda::ccmetagen=1.1.4"
    // container "${ workflow.containerEngine == 'singularity' ?
    //     'oras://community.wave.seqera.io/library/ccmetagen:1.1.4' :
    //     'biocontainers/ccmetagen:1.1.4' }"

    input:
    tuple val(meta), path(res), path(mapstat)
    tuple val(db_id), path(db_dir)

    output:
    tuple val(meta), path("${meta.id}_ccmetagen/")      , emit: results_dir
    tuple val(meta), path("${meta.id}_ccmetagen/*.csv") , emit: csv
    path  "versions.yml"                                , emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Locate the DB prefix in the staged DB directory
    // CCMetagen needs it for NCBI taxonomy lookups
    """
    DB_PREFIX=\$(find "${db_dir}" \\
        -maxdepth 2 -name "*.name" 2>/dev/null \\
        | head -1 | sed 's/\\.name\$//')

    if [[ -z "\$DB_PREFIX" ]]; then
        echo "[WARN] Could not locate KMA DB prefix — running without --db"
        DB_FLAG=""
    else
        DB_FLAG="--db \$DB_PREFIX"
    fi

    mkdir -p ${prefix}_ccmetagen

    CCMetagen.py \\
        -i "${res}" \\
        -o ${prefix}_ccmetagen/${prefix} \\
        \$DB_FLAG \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CCMetagen: \$(CCMetagen.py --version 2>&1 | head -1 | sed 's/CCMetagen //' || echo "N/A")
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_ccmetagen
    echo "sample,kingdom,phylum,class,order,family,genus,species,reads" \\
        > ${meta.id}_ccmetagen/${meta.id}.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CCMetagen: 1.1.4
        python: 3.11.0
    END_VERSIONS
    """
}

process KMA {
    tag "$sample_id"

    // container "ccmetagen_kma.sif"

    input:
    tuple val(sample_id), path(reads)
    tuple val(db_id), path(db_dir)

    output:
    tuple val(sample_id), path("${sample_id}.res"), path("${sample_id}.mapstat")

    script:
    """
    kma \
        -ipe ${reads[0]} ${reads[1]} \
        -o ${sample_id} \
        -t_db ${db_dir} \
        -t $task.cpus

    mv ${sample_id}.res ${sample_id}.res
    mv ${sample_id}.mapstat ${sample_id}.mapstat
    """
}


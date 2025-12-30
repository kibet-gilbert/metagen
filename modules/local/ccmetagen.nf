process CCMETAGEN {
    tag "$sample_id"

    // container "ccmetagen_ccmetagen.sif"

    input:
    tuple val(sample_id), path(res), path(mapstat)
    tuple val(db_id), path(db_dir)

    output:
    tuple val(sample_id), path("${sample_id}_ccmetagen")

    script:
    """
    mkdir ${sample_id}_ccmetagen

    CCMetagen.py \
        -i ${res} \
        -o ${sample_id}_ccmetagen \
        --db ${db_dir}

    """
}


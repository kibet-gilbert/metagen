
process RESFINDER {
    tag "${sample_id} + ${ResFinderName}"

    input:
    tuple val(sample_id), path(nohost_reads)                         
    tuple val(ResFinderName), path(RESFINDER_db)
    tuple val(PointFinderName), path(POINTFINDER_db)
    tuple val(DisinFinderName), path(DISINFINDER_db)

    output:
    tuple val(sample_id), path("*.txt"), emit: results
    tuple val(sample_id), path("*.json"), emit: results_json
    tuple val(sample_id), path("*_seq.fsa"), emit: sequences
    tuple val(sample_id), path("${sample_id}/"), emit: resDir

    script:
    """
    # Setting Databases environment variables
    CGE_RESFINDER_RESGENE_DB=${RESFINDER_db}
    CGE_RESFINDER_RESPOINT_DB=${POINTFINDER_db}
    CGE_DISINFINDER_DB=${DISINFINDER_db}
    export MPLCONFIGDIR=\$(pwd)/.cache/matplotlib
    export TMPDIR=\$(pwd)/tmp
    mkdir -p "\$TMPDIR"
    mkdir -p "\$MPLCONFIGDIR"

    # Clean any residual *_kma folders:
    rm -rf *_kma/

    python3 -m resfinder \\
        --inputfastq ${nohost_reads[0]} ${nohost_reads[1]} \\
        --species "other" \\
        --min_cov 0.6 \\
        --threshold 0.9 \\
        --acq_overlap 30 \\
        --outputPath ./ \\
        --out_json ./${sample_id}.resfinder.json \\
        --acquired \\
        --db_path_res ${RESFINDER_db} \\
        --point \\
        --min_cov_point 0.6 \\
        --threshold_point 0.9 \\
        --db_path_point ${POINTFINDER_db} \\
        --disinfectant \\
        --db_path_disinf ${DISINFINDER_db} \\
        --ignore_stop_codons \\
        --ignore_indels \\
        --kma_threads $task.cpus 

    # Rename the output files to include ${sample_id}
    for file in *.txt *.json *_seq.fsa; do
        if [[ -e "\$file" ]]; then
            mv "\$file" "${sample_id}.\$file"
        fi
    done
    mkdir -p ${sample_id}
    cp -r disinfinder_kma/ ${sample_id}/
    cp -r pointfinder_kma/ ${sample_id}/
    cp -r resfinder_kma/ ${sample_id}/
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       resfinder: \$(echo \$( python -m resfinder -v 2>&1) )
    END_VERSIONS
    """
}

process CENTRIFUGE {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    tuple val(centrifugeDBName), path("database/*")

    output:
    tuple val(sample_id), path("*_report.txt"), emit: centrifuge_report
    tuple val(sample_id), path("*_results.txt"), emit: centrifuge_results
    tuple val(sample_id), path("*_kreport2.txt"), emit: kraken_report
    tuple val(sample_id), path("*_centrifuge.krona"), emit: krona_report

    script:
    centrifugeSIFRun='apptainer run /home/gkibet/bioinformatics/github/metagenomics/scripts/apptainer/centrifuge_imp_python.sif'
    """
    #${centrifugeSIFRun} \\
    centrifuge -x database/"${centrifugeDBName}" \\
           -p ${task.cpus} \\
           --report-file ${sample_id}_report.txt \\
           -S ${sample_id}_results.txt \\
           -1 ${trimmed_reads[0]} \\
           -2 ${trimmed_reads[1]}
           #-r ${trimmed_reads[2]}
    #${centrifugeSIFRun} \\
    centrifuge-kreport -x database/"${centrifugeDBName}" \\
           ${sample_id}_results.txt > ${sample_id}.centrifuge_kreport2.txt \\
           2> >(tee ${sample_id}.centrifuge_kreport2_skipped.txt >&2)
    cat ${sample_id}_results.txt | cut -f 1,3 > ${sample_id}_centrifuge.krona
    """
}

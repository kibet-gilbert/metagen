process FASTP {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_trim*"), emit: trimmed_reads
    tuple val(sample_id), path("*failed*"), emit: failed_reads
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.json"), emit: json
    tuple val(sample_id), path("*.log"), emit: log

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${sample_id}_R1_trim.fastq.gz \\
        --out2 ${sample_id}_R2_trim.fastq.gz \\
        --html ${sample_id}.fastp.html \\
	--json ${sample_id}.fastp.json \\
        --failed_out ${sample_id}_failed.fastq.gz \\
        --thread $task.cpus \\
        --detect_adapter_for_pe \\
        --dedup \\
        |& tee -a -i ${sample_id}.fastp.log

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
          fastp: \$(echo \$(fastp -v 2>&1) | sed 's/fastp //')
    END_VERSIONS
    """
}

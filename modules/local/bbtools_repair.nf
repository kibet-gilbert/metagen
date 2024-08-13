process BBTOOLSREPAIR {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(samFile)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: trimmed_reads

    script:
    """
    repair.sh in=${samFile} \\
        out1=${sample_id}.argSeqs.R1.fastq \\
        out2=${sample_id}.argSeqs.R2.fastq \\
        outs=${sample_id}.argSeqs.SE.fastq \\
        repair=t
    pigz ${sample_id}.argSeqs.R1.fastq
    pigz ${sample_id}.argSeqs.R2.fastq
    pigz ${sample_id}.argSeqs.SE.fastq
    """
}


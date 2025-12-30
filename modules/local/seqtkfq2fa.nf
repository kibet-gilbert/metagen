process FASTQ2FASTA {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)  // expects [R1.fastq.gz, R2.fastq.gz]

    output:
    tuple val(sample_id), path("${sample_id}_R1.fasta"), emit: fasta_r1
    tuple val(sample_id), path("${sample_id}_R2.fasta"), emit: fasta_r2

    script:
    """
    echo "Converting FASTQ to FASTA for: ${sample_id}"

    seqtk seq -a ${reads[0]} > ${sample_id}_R1.fasta
    seqtk seq -a ${reads[1]} > ${sample_id}_R2.fasta

    """
}

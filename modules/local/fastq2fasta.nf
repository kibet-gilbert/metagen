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

    gunzip -c ${reads[0]} | \\
        awk 'NR%4==1 {printf(">%s\\n",substr(\$0,2))} NR%4==2 {print}' \\
        > ${sample_id}_R1.fasta

    gunzip -c ${reads[1]} | \\
        awk 'NR%4==1 {printf(">%s\\n",substr(\$0,2))} NR%4==2 {print}' \\
        > ${sample_id}_R2.fasta
    """
}

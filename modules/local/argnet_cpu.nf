process ARGNET_CPU {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(input_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.argnet.tsv"), emit: argnet_output

    script:
    """
        python -m argnet \\
        --input ${input_fasta} \\
        --type nt \\
        --model 'argnet-s' \\
        --outname ${sample_id}.argnet.tsv
    """
}


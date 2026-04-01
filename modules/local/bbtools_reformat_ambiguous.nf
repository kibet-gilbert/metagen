process BBTOOLSREFORMAT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(nohost_reads)    // [R1, R2]

    output:
    tuple val(sample_id), path("*.noambig.*"), emit: maskedreads

    script:
    def r1 = nohost_reads[0]
    def r2 = nohost_reads[1]

    """
    reformat.sh in=${r1} in2=${r2} \\
        out1=${sample_id}.noambig.R1.fastq.gz \\
        out2=${sample_id}.noambig.R2.fastq.gz \\
        iupacToN=t \\
        tuc=t \\
        fixjunk=t \\
        dotdashxton=t
    """
}


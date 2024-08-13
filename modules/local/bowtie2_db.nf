/*
 * Bowtie2 for index build
 */
process BOWTIE2_DB {
    tag "$PhiXAccession"

    input:
    tuple val(PhiXAccession), path(fasta)

    output:
    tuple val(PhiXAccession), path('bt2_index_base*'), emit: index

    script:
    """
    bowtie2-build --threads ${task.cpus} ${fasta} bt2_index_base
    """
}

process KRAKEN2FILTER {
    // publishDir "${params.outputDir}/kraken2", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(trimmed_reads)
    path("database/*")

    output:
    tuple val(sample_id), path("*.nohost{.,_}*"), emit: nohost_reads
    tuple val(sample_id), path("*.report.txt"), emit: kraken_hostmap_report1
    tuple val(sample_id), path("*.kraken2.out"), emit: kraken_hostmap_report2

    script:
    """
    kraken2 --db database/ --threads $task.cpus \\
            --unclassified-out ${sample_id}.nohost#.fastq \\
            --classified-out ${sample_id}.host#.fastq \\
            --report ${sample_id}.kraken2.report.txt \\
            --output ${sample_id}.kraken2.out \\
            --gzip-compressed --report-zero-counts \\
            --paired $trimmed_reads
    """
}


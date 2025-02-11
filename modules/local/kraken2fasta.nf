process KRAKEN2FASTA {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(contigs)
    path("database/*")
    // path(kraken2_db)

    output:
    tuple val(sample_id), path("*_kreport.txt"), emit: kraken_kreport
    tuple val(sample_id), path("*_kraken2.out"), emit: kraken_results
    tuple val(sample_id), path("*_kraken2.krona"), emit: krona_report

    script:
    """
    kraken2 --db database/ \\
            --threads $task.cpus \\
            --report ${sample_id}_taxContigs_kreport.txt \\
            --output ${sample_id}_taxContigs_kraken2.out \\
            --report-zero-counts \\
            ${contigs}
    cat ${sample_id}_taxContigs_kraken2.out | cut -f 2,3 > ${sample_id}_taxContigs_kraken2.krona
    """
}


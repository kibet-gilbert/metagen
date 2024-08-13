process KRAKEN2 {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    path("database/*")

    output:
    //tuple val(sample_id), path("*.classified{.,_}*"), emit: classified_reads
    tuple val(sample_id), path("*_kreport.txt"), emit: kraken_kreport
    tuple val(sample_id), path("*_kraken2.out"), emit: kraken_report
    tuple val(sample_id), path("*_kraken2.krona"), emit: krona_report

    script:
    """
     kraken2 --db database/ \\
            --threads $task.cpus \\
            --unclassified-out ${sample_id}.unclassified#.fastq \\
            --classified-out ${sample_id}.classified#.fastq \\
            --report ${sample_id}_tax_kreport.txt \\
            --output ${sample_id}_tax_kraken2.out \\
            --report-zero-counts \\
            --gzip-compressed \\
            --paired ${trimmed_reads[0]} ${trimmed_reads[1]}
    cat ${sample_id}_tax_kraken2.out | cut -f 2,3 > ${sample_id}_kraken2.krona
    """
}

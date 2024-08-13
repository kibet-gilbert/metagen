process KRONA {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(krona_report)                                   
    path(taxonomy)

    output:
    tuple val(sample_id), path("*.html"), emit: krona_html

    script:
    """
    ktImportTaxonomy -tax ${taxonomy} \\
                     -o ${sample_id}_Kraken2.krona.html \\
                     ${krona_report}
    """
}


process CARDRGI {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(trimmed_reads)                                  
    path(RGI_db)

    output:
    tuple val(sample_id), path("*allele_mapping_data.txt"), emit: allele
    tuple val(sample_id), path("*gene_mapping_data.txt"), emit: gene
    tuple val(sample_id), path("*.temp.sam.temp.*"), emit: 
    tuple val(sample_id), path("*.sorted.length_100.*"), emit: 
    tuple val(sample_id), path("*overall_mapping_stats*"), emit: 
    tuple val(sample_id), path("*.sorted.length_100*"), emit: 

    script:
    """
    rgi bwt \\
        -1 ${trimmed_reads[0]} \\
        -2 ${trimmed_reads[1]} \\
        -a kma \\
        -n $task.cpus \\
        -o ${sample_id} \\
        --local ${RGI_db} \\
        --include_other_models \\
        --include_wildcard
    """
}


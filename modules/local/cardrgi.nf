
process CARDRGI {
    tag "${sample_id}+${RGIDBName}"

    input:
    tuple val(sample_id), path(trimmed_reads)
    tuple val(RGIDBName), path(rgiDBpath) 
    // path(RGI_db)

    output:
    tuple val(sample_id), path("*allele_mapping_data.txt"), emit: allele
    tuple val(sample_id), path("*gene_mapping_data.txt"), emit: gene
    tuple val(sample_id), path("*.sorted.length_100.bam"), emit: argbamfile
    tuple val(sample_id), path("*.sorted.length_100.bam.bai"), emit: argbamindex
    tuple val(sample_id), path("*overall_mapping_stats*"), emit: mapping_stats

    script:
    """
    rgi bwt \\
        -1 ${trimmed_reads[0]} \\
        -2 ${trimmed_reads[1]} \\
        -a kma \\
        -n $task.cpus \\
        -o ${sample_id} \\
        --local \\
        --include_other_models \\
        --include_wildcard
    """
}

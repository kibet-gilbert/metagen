process CARDRGI {
    tag "${sample_id}+${RGIDBName}"

    input:
    tuple val(sample_id), path(trimmed_reads)
    tuple val(RGIDBName), path(rgiDBpath) 
    // path(RGI_db)

    output:
    tuple val(sample_id), path("*allele_mapping_data.txt"), emit: allele
    tuple val(sample_id), path("*gene_mapping_data.txt"), emit: gene
    tuple val(sample_id), path("*.allele_mapping_data.json"), emit: argjsonfile
    // tuple val(sample_id), path("*_mapping_stats.txt"), emit: argbamindex
    tuple val(sample_id), path("*_mapping_stats.txt"), emit: mapping_stats

    script:
    """
    # Override HOME and MPLCONFIGDIR to point to this local work directory
    export HOME=\$(pwd)/.fake_home
    export MPLCONFIGDIR=\$HOME/.matplotlib
    export XDG_CACHE_HOME=\$HOME/.cache
    export TMPDIR=\$(pwd)/.tmp
    mkdir -p \$MPLCONFIGDIR  \$TMPDIR \$XDG_CACHE_HOME

    # Run RGI
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

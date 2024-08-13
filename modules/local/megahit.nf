process MEGAHIT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(nohost_reads)                                   

    output:
    tuple val(sample_id), path("${sample_id}.contigs.fa"), emit: contigs       
    tuple val(sample_id), path("*.log"), emit: megahitlog
    tuple val(sample_id), path("${sample_id}.contigs.fa.gz"), emit: assembly_gz
    path("versions.yml"), emit: versions

    script:
    mem = task.memory.toBytes()
    """
    megahit \\
        -t "${task.cpus}" \\
        -m $mem \\
        -1 ${nohost_reads[0]} \\
        -2 ${nohost_reads[1]} \\
        -o . \\
        --out-prefix "${sample_id}"
    gzip -c "${sample_id}.contigs.fa" > "${sample_id}.contigs.fa.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """
}


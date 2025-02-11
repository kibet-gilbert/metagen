/*
 * hostile for read removal
 */
process HOSTILE {
    tag "${sample_id}+${hostileDBName}"

    // conda "bioconda::hostile=1.1.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ef/ef586e13547ecd7533cef8ae7b663cbf7425dd13b52180481b882e2e10e7f544/data' :
    //     'oras://community.wave.seqera.io/library/hostile:1.1.0--f22a8f42b4582307' }"

    input:
    tuple val(sample_id), path(reads)
    tuple val(hostileDBName), path(hostileDBpath)

    output:
    tuple val(sample_id), path("*.clean*.fastq.gz"), emit: reads
    // path "*.mapped*.read_ids.txt", optional:true, emit: read_ids
    tuple val(sample_id), path("*.hostile.logs.json"), emit: log

    script:
    // println "hostileDBName: ${hostileDBName}"
    hostileSIFRun='apptainer run /home/gkibet/bioinformatics/github/apptainer/hostile/hostile_apptainer.sif'
    """
    export HOSTILE_CACHE_DIR=./${hostileDBpath}/
    #${hostileSIFRun} \\
    hostile clean \\
        --offline \\
        --out-dir . \\
    	--threads ${task.cpus} \\
   	--fastq1 "${reads[0]}" \\
   	--fastq2 "${reads[1]}" \\
   	--index ${hostileDBName} \\
        |& tee -a -i ${sample_id}.hostile.logs.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       hostile: \$(echo \$( hostile --version 2>&1) )
    END_VERSIONS
    """
}

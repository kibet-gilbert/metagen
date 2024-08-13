
process FASTQC {
    // publishDir "${params.outputDir}/fastqc", mode: 'copy'
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip

    script:
    """
    fastqc -t $task.cpus -o . ${reads}
    """
}



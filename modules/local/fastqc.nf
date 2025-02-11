
process FASTQC {
    // publishDir "${params.outputDir}/fastqc", mode: 'copy'
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip

    // container '/home/gkibet/bioinformatics/github/apptainer/fastqc/fastqc_apptainer.sif'

    script:
    fastqcSIFRun='apptainer run /home/gkibet/bioinformatics/github/apptainer/fastqc/fastqc_apptainer.sif'
    """
    #${fastqcSIFRun} \\
    fastqc -t $task.cpus -o . ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       fastqc: \$(echo \$( fastqc --version 2>&1) | sed 's/FastQC v//' )
    END_VERSIONS
    """
}



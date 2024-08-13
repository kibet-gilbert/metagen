process SAM2FASTQ {
    // debug true
    tag "$sample_id"

    input:
    tuple val(sample_id), path(argbamfile)                                     

    output:
    path("*.flagstat.txt"), emit: flagstat_report                              
    tuple val(sample_id), path("*.fastq.gz"), emit: reads                      
    tuple val(sample_id), path("*.sam"), emit: samFile                         
    tuple val(sample_id), path("*.log"), emit: log                             

    script:
    //myBamFile = file()
    """
    samtools flagstat ${argbamfile} \\
        -@ $task.cpus \\
        |& tee ${sample_id}.flagstat.txt
    samtools view ${argbamfile} \\
        -hu \\
        -G 0xd \\
        -@ $task.cpus |
    samtools sort -n - \\
        -@ $task.cpus \\
        -o ${sample_id}.sortednFiltered.argSeqs.sam \\
        |& tee -a ${sample_id}.samtoolsfastq.log
    samtools fastq ${sample_id}.sortednFiltered.argSeqs.sam \\
        -@ $task.cpus \\
        -1 ${sample_id}.argSeqs_R1.fastq.gz \\
        -2 ${sample_id}.argSeqs_R2.fastq.gz \\
        -n \\
        |& tee -a ${sample_id}.samtoolsfastq.log
        ##-0 ${sample_id}.argSeqs_RO.fastq.gz \\
        ##-s ${sample_id}.argSeqs_SR.fastq.gz \\
    """
}


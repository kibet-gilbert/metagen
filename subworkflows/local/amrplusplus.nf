nextflow.enable.dsl=2

workflow AMRPLUSPLUS {
    take:
    tuple val(sample_id), path(reads)

    main:
    call AMRplusplus with: sample_id, reads

    emit:
    amr_output = AMRplusplus.out
}

process AMRplusplus {
    tag "$sample_id"

    input:
    val sample_id
    path reads

    output:
    path "amrplusplus_${sample_id}", emit: out

    container:
    'apptainer/amrplusplus.sif'

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    mgs2amr.sh \\
        -i ${read1} \\
        -j ${read2} \\
        -o ./ \\
        -n ${sample_id} \\
        -c $task.cpus \\
        -m 32G \\
        -v 1
    """
}


process MGS2AMR {

    tag "${sample_id}+mgs2amr_db+${blastdb_name}"
    containerOptions "--bind ${params.blast_db}:${params.blast_db},\$(pwd):/opt/mgs2amr/dataAndScripts"

    input:
    tuple val(sample_id), path(reads)        // [R1, R2]
    path dataAndScripts_dir
    tuple val(blastdb_name), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}/"), emit: amr_results

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    def step_flag   = params.mgs2amr_step    ? "-s ${params.mgs2amr_step}" : ""
    def verbose_flag = params.mgs2amr_verbose ? "-v ${params.mgs2amr_verbose}" : ""
    def force_flag  = params.mgs2amr_force   ? "-f" : ""
    def compress_flag = params.mgs2amr_compress ? "-z TRUE" : "-z FALSE"
    // def db_flag      = db_path  ? "-d ${db_path}"    : ""
    def mem_flag = "-m ${task.memory.toGiga()}G"

    """
    export TMPDIR=\$(pwd)/tmp
    export HOME=\$(pwd)/.fake_home
    export BLASTDB=${blast_db}
    export DATAANDSCRIPTS=\$(pwd)
    export _JAVA_OPTIONS="-XX:CompressedClassSpaceSize=256m -XX:MaxMetaspaceSize=512m"
    mkdir -p \${TMPDIR} \${HOME} #\${DATAANDSCRIPTS}
    for f in ${dataAndScripts_dir}/*; do
        cp -r "\$(readlink -f "\$f")" "\${DATAANDSCRIPTS}/"
        # chmod 777 \${DATAANDSCRIPTS}/*
    done

    mgs2amr.sh \\
      -i ${r1} \\
      -j ${r2} \\
      -o ./ \\
      -n ${sample_id} \\
      -c $task.cpus \\
      -d "./mgs2amr.db" \\
      ${mem_flag} \\
      ${step_flag} \\
      ${verbose_flag} \\
      ${force_flag} \\
      ${compress_flag}
    """
}


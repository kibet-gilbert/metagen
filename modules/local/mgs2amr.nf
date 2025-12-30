process MGS2AMR {

    tag "${sample_id}" // +${db_name}"
    containerOptions "--bind \$(pwd):/opt/mgs2amr/dataAndScripts"

    input:
    tuple val(sample_id), path(reads)        // [R1, R2]
    val step                                 // e.g. 4
    val verbose                              // e.g. 1
    val compress                             // e.g. true
    val force                                // e.g. true
    tuple val(db_name), path(db_path)

    output:
    tuple val(sample_id), path(".*tar.gz"), emit: amr_results

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    def step_flag    = step     ? "-s ${step}"       : ""
    def verbose_flag = verbose  ? "-v ${verbose}"    : ""
    def force_flag   = force    ? "-f"               : ""
    def compress_flag= compress ? "-z TRUE"          : "-z FALSE"
    def db_flag      = db_path  ? "-d ${db_path}"    : ""

    """
    export TMPDIR=\$(pwd)/tmp
    export HOME=\$(pwd)/.fake_home
    export BLASTDB=${params.blast_db}
    mkdir -p \$TMPDIR \$HOME

    mgs2amr.sh \\
      -i ${r1} \\
      -j ${r2} \\
      -o ./ \\
      -n ${sample_id} \\
      -c $task.cpus \\
      -m $task.memory \\
      ${db_flag}
      ${step_flag} \\
      ${verbose_flag} \\
      ${force_flag} \\
      ${compress_flag}
      # ${db_flag}
      # -d ./mgs2amr.db \\
    """
}


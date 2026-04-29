// ============================================================================
// modules/local/kma_align.nf
//
// Align paired-end reads against a KMA-indexed database and output:
//   - .res  : per-template summary (identity, coverage, depth, reads mapped)
//   - .frag : per-fragment (read-pair) alignment details (for abundance)
//   - .aln  : consensus alignment
//   - .fsa  : consensus FASTA
//   - .mat  : per-position depth matrix (optional, large)
//
// This module is parameterised for three distinct use cases:
//   mode = 'scg'       : GTDB bac120 SCGs — GEQ normalisation
//                        relaxed identity (80%), includes fragment output
//   mode = 'ccmetagen' : CCMetagen classification (RefSeq / core_nt / nt)
//                        standard identity, CCMetagen-compatible output
//   mode = 'resfinder' : ARG quantification (ResFinder/CARD)
//                        strict identity (90%), matches existing pipeline
//
// The mode is set via task.ext.mode or params.kma_mode.
// All modes produce a .res file in identical column format for easy merging.
//
// Input channel:
//   tuple val(meta), path(reads), val(db_name), path(kma_index_dir)
//     meta         : map with keys: id, sample, county, site, week, total_reads
//     reads        : list [R1.fastq.gz, R2.fastq.gz]
//     db_name      : database label (e.g. 'bac120_scg_r232')
//     kma_index_dir: directory containing KMA index files
//
// Output channels:
//   res    : tuple val(meta), path("*.res")         — primary abundance table
//   frag   : tuple val(meta), path("*.frag.gz")     — per-read alignment
//   aln    : tuple val(meta), path("*.aln")         — consensus alignment
//   fsa    : tuple val(meta), path("*.fsa")         — consensus FASTA
//   log    : tuple val(meta), path("*.kma.log")     — run log
//
// Usage:
//   include { KMA_ALIGN } from './modules/local/kma_align'
//
//   // SCG / GEQ mode:
//   KMA_ALIGN(HOSTILE.out.reads,kma_index)
//
//   // CCMetagen mode:
//   KMA_ALIGN(HOSTILE.out.reads,ccmetagenDB_index)
// ============================================================================

process KMA_ALIGN {
    tag "${meta.id}-${db_name}"
    label 'process_medium'

    // ── Container ─────────────────────────────────────────────────────────
    conda     "bioconda::kma=1.6.8"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/kma:1.6.8--h577a1d6_0' :
        'biocontainers/kma:1.6.8--h577a1d6_0' }"

    // ── Publish ───────────────────────────────────────────────────────────
    publishDir(
        path:    "${params.outdir}/kma/${meta.id}/${db_name}",
        mode:    params.publish_dir_mode ?: 'copy',
        pattern: "*.{res,frag.gz,aln,fsa,kma.log}"
    )

    // ── Inputs ────────────────────────────────────────────────────────────
    input:
    tuple val(meta), path(reads)
    tuple val(db_name), path(kma_index_dir)

    // ── Outputs ───────────────────────────────────────────────────────────
    output:
    tuple val(meta), path("${prefix}.res"),      emit: res
    tuple val(meta), path("${prefix}.frag.gz"),  emit: frag,    optional: true
    tuple val(meta), path("${prefix}.aln"),      emit: aln,     optional: true
    tuple val(meta), path("${prefix}.fsa"),      emit: fsa,     optional: true
    tuple val(meta), path("${prefix}.kma.log"),  emit: log
    path  "versions.yml",                        emit: versions

    // ── Script ────────────────────────────────────────────────────────────
    script:
    prefix         = "${meta.id}__${db_name}"
    def db_prefix  = "${kma_index_dir}/${db_name}"
    def r1         = reads[0]
    def r2         = reads[1]
    def mode       = task.ext.mode ?: params.kma_mode ?: 'scg'
    def extra_args = task.ext.args ?: ''

    // ── Mode-specific flag sets ────────────────────────────────────────────
    // scg mode: relaxed identity for divergent marker genes
    //           fragment output enabled for per-read abundance
    //           lower ID threshold captures reads from organisms divergent
    //           from the representative genome
    //    -ID 80 \\
    //    -coverage 80 \\
    def scg_flags = """\
        -1t1 \\
        -mem_mode \\
        -and \\
        -apm f \\
        -ef \\
        -ID ${params.kma_identity ?: 80} \\
        -coverage ${params.kma_coverage ?: 80} \\
        -t ${task.cpus}"""

    // ccmetagen mode: standard KMA settings for CCMetagen pipeline
    //                 -sam disabled (too large for metagenomics)
    //                 -ef enables fragment output for abundance
    //    -ID 50 \\
    //    -coverage 50 \\
    def ccmetagen_flags = """\
        -1t1 \\
        -mem_mode \\
        -and \\
        -apm f \\
        -ef \\
        -ID ${params.kma_identity ?: 50} \\
        -coverage ${params.kma_coverage ?: 50} \\
        -t ${task.cpus}"""

    // resfinder mode: strict identity matching ARG pipeline settings
    //                 matches existing ResFinder/CARD-RGI workflow
    //    -ID 90 \\
    //    -coverage 90 \\
    def resfinder_flags = """\
        -1t1 \\
        -mem_mode \\
        -and \\
        -apm f \\
        -ef \\
        -ID ${params.kma_identity ?: 90} \\
        -coverage ${params.kma_coverage ?: 90} \\
        -t ${task.cpus}"""

    def mode_flags = (mode == 'ccmetagen') ? ccmetagen_flags :
                     (mode == 'resfinder') ? resfinder_flags :
                     scg_flags   // default to scg

    """
    echo "[INFO] KMA alignment: ${meta.id} vs ${db_name} (mode: ${mode})"
    echo "[INFO] R1: ${r1}"
    echo "[INFO] R2: ${r2}"
    echo "[INFO] DB: ${db_prefix}"

    # ── Run KMA ─────────────────────────────────────────────────────────
    kma \\
        -ipe ${r1} ${r2} \\
        -o ${prefix} \\
        -t_db ${db_prefix} \\
        ${mode_flags} \\
        ${extra_args} \\
        2>&1 | tee ${prefix}.kma.log

    # ── Validate output ──────────────────────────────────────────────────
    if [[ ! -f "${prefix}.res" ]]; then
        echo "[ERROR] KMA produced no .res file — alignment may have failed" >&2
        exit 1
    fi

    n_templates=\$(tail -n +2 ${prefix}.res | wc -l)
    echo "[INFO] Templates with reads mapped: \${n_templates}"

    # ── Compress fragment file if present ────────────────────────────────
    if [[ -f "${prefix}.frag" ]]; then
        gzip -f ${prefix}.frag
        echo "[INFO] Fragment file compressed: ${prefix}.frag.gz"
    fi

    # ── Append sample metadata to .res header for downstream merging ─────
    # Inserts a comment line with meta fields before the column header
    meta_line="# sample=${meta.id} db=${db_name} mode=${mode} total_reads=${meta.total_reads ?: 'NA'} county=${meta.county ?: 'NA'} site=${meta.site ?: 'NA'} week=${meta.week ?: 'NA'}"
    sed -i "1i\\\${meta_line}" ${prefix}.res

    echo "[INFO] Alignment complete: ${prefix}.res"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(kma -v 2>&1 | grep -oP 'KMA-\\K[0-9.]+' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    prefix = "${meta.id}__${db_name}"
    """
    # Minimal stub .res matching KMA column format
    cat <<-EOF > ${prefix}.res
    # sample=${meta.id} db=${db_name} mode=stub
    #Template\\tScore\\tExpected\\tTemplate_length\\tTemplate_Identity\\tTemplate_Coverage\\tQuery_Identity\\tQuery_Coverage\\tDepth\\tq_value\\tp_value
    stub_template|stub_marker\\t100\\t95\\t300\\t99.0\\t100.0\\t99.0\\t100.0\\t5.2\\t0.0\\t0.0
    EOF
    touch ${prefix}.frag.gz
    touch ${prefix}.aln
    touch ${prefix}.fsa
    touch ${prefix}.kma.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: 1.6.8
    END_VERSIONS
    """
}

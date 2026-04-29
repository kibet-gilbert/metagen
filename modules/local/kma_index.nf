// ============================================================================
// modules/local/kma_index.nf
//
// Build a KMA index from a FASTA file.
// Supports any FASTA source:
//   - GTDB bac120 SCG (from build_gtdb_scg_kma.sh)
//   - CCMetagen RefSeq / core_nt / nt (from download_all_refseq.sh or clean_fasta.sh)
//   - Any other nucleotide FASTA
//
// Input channel:
//   tuple val(db_name), path(fasta)
//     db_name : short label used for output directory and file naming
//     fasta   : compressed (.gz) or plain FASTA to be indexed
//
// Output channel:
//   tuple val(db_name), path("${db_name}_kma_index")
//     db_name         : passed through for downstream use
//     kma_index_dir   : directory containing all KMA index files
//
// Usage in workflow:
//   include { KMA_INDEX } from './modules/local/kma_index'
//
//   Channel
//     .of(
//       [ 'bac120_scg_r232',  file(params.gtdb_scg_fasta)   ],
//       [ 'refseq_ccmetagen', file(params.refseq_fasta)      ],
//       [ 'core_nt_ccmetagen',file(params.core_nt_fasta)     ]
//     )
//     .set { ch_fastas_to_index }
//
//   KMA_INDEX(ch_fastas_to_index)
// ============================================================================

process KMA_INDEX {
    tag "${db_name}"
    label 'process_high_memory'

    // ── Container ─────────────────────────────────────────────────────────
    conda     "bioconda::kma=1.6.8"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/kma:1.6.8--h577a1d6_0' :
        'biocontainers/kma:1.6.8--h577a1d6_0' }"

    // ── Publish ───────────────────────────────────────────────────────────
    publishDir(
        path:    "${params.outdir}/kma_databases/${db_name}",
        mode:    params.publish_dir_mode ?: 'copy',
        pattern: "${db_name}_kma_index/**"
    )

    // ── Inputs ────────────────────────────────────────────────────────────
    input:
    tuple val(db_name), path(fasta)

    // ── Outputs ───────────────────────────────────────────────────────────
    output:
    tuple val(db_name), path("${db_name}_kma_index"), emit: index
    path  "versions.yml",                             emit: versions

    // ── Script ────────────────────────────────────────────────────────────
    script:
    def args       = task.ext.args ?: ''
    def prefix     = "${db_name}_kma_index/${db_name}"
    def input_flag = fasta.name.endsWith('.gz') ? "-i ${fasta}" : "-i ${fasta}"
    // KMA index handles gzipped FASTA natively

    """
    mkdir -p ${db_name}_kma_index

    echo "[INFO] Building KMA index: ${db_name}"
    echo "[INFO] Input FASTA: ${fasta}"
    echo "[INFO] Index prefix: ${prefix}"

    # Count input sequences for logging
    if [[ "${fasta}" == *.gz ]]; then
        n_seq=\$(zgrep -c '^>' "${fasta}" || echo "unknown")
    else
        n_seq=\$(grep  -c '^>' "${fasta}" || echo "unknown")
    fi
    echo "[INFO] Input sequences: \${n_seq}"

    # Build the KMA index
    kma index \\
        ${input_flag} \\
        -o ${prefix} \\
        ${args} \\
        2>&1 | tee ${db_name}_kma_index/kma_index.log

    # Verify index files were created
    expected_files="${prefix}.index ${prefix}.seq ${prefix}.name ${prefix}.length.b ${prefix}.comp.b"
    for f in \${expected_files}; do
        if [[ ! -f "\${f}" ]]; then
            echo "[ERROR] Expected index file missing: \${f}" >&2
            exit 1
        fi
    done

    echo "[INFO] KMA index built successfully"
    echo "[INFO] Index files:"
    ls -lh ${db_name}_kma_index/ | tee -a ${db_name}_kma_index/kma_index.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(kma -v 2>&1 | grep -oP 'KMA-\\K[0-9.]+' || kma 2>&1 | head -1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${db_name}_kma_index
    touch ${db_name}_kma_index/${db_name}.index
    touch ${db_name}_kma_index/${db_name}.seq
    touch ${db_name}_kma_index/${db_name}.name
    touch ${db_name}_kma_index/${db_name}.length.b
    touch ${db_name}_kma_index/${db_name}.comp.b
    touch ${db_name}_kma_index/${db_name}.decon.comp.b

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: 1.6.8
    END_VERSIONS
    """
}

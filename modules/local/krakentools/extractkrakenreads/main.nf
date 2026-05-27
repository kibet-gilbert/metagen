// =============================================================================
// modules/local/krakentools/extractkrakenreads/main.nf
//
// Extracts reads from FASTQ files based on Kraken2/Bracken classification
// using KrakenTools' extract_kraken_reads.py.
//
// Use case in pathogen validation workflow:
//   Given Kraken2 classifications, pull out the reads assigned to a list of
//   target pathogen taxIDs (and optionally their taxonomic children), then
//   pass those reads to BLAST/MAG assembly for orthogonal confirmation.
//
// Inputs:
//   tuple val(meta), path(reads), path(kraken_output), path(kraken_report)
//     meta            : [id: <sample_id>, ...] standard nf-core meta map
//     reads           : [R1.fastq.gz, R2.fastq.gz] OR [R1.fastq.gz] for SE
//     kraken_output   : raw per-read kraken2 output file (NOT the report)
//     kraken_report   : kraken2 .report file (required for --include-children)
//   val(taxids)       : list of taxIDs (integers or strings) to extract
//
// Outputs:
//   tuple val(meta), path("*.extracted_{1,2}.fasta.gz") : extracted reads
//   tuple val(meta), path("*.summary.txt")              : extraction summary
//   path "versions.yml"
// =============================================================================

process KRAKENTOOLS_EXTRACTKRAKENREADS {
    tag "${meta.id}"
    label 'process_medium'

    conda     "bioconda::krakentools=1.2 conda-forge::pigz=2.8"
    container "${ (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer')?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0' :
        'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0' }"
        // 'https://depot.galaxyproject.org/singularity/mulled-v2-ef891e9b0617d884630ceb63401a526b6008baec:5f7c6be2f66db88a16bfdb900e71a9800ddbcb04' :
        // 'quay.io/biocontainers/mulled-v2-ef891e9b0617d884630ceb63401a526b6008baec:5f7c6be2f66db88a16bfdb900e71a9800ddbcb04' }"

    publishDir(
        path:  { "${params.outdir}/krakentools/extract_reads/${meta.id}" },
        mode:  params.publish_dir_mode ?: 'copy',
        pattern: "*.{fasta.gz,fastq.gz,summary.txt}"
    )

    input:
    tuple val(meta), path(reads), path(kraken_output), path(kraken_report)
    val   taxids

    output:
    tuple val(meta), path("${prefix}.extracted_1.${suffix}.gz"),                    emit: extracted_r1
    tuple val(meta), path("${prefix}.extracted_2.${suffix}.gz"),  optional: true,  emit: extracted_r2
    tuple val(meta), path("${prefix}.summary.txt"),                                 emit: summary
    path  "versions.yml",                                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ── Resolve parameters ────────────────────────────────────────────────
    prefix         = task.ext.prefix ?: "${meta.id}"
    def args       = task.ext.args   ?: ''
    def fastq_out  = task.ext.fastq_output != null ? task.ext.fastq_output : false
    suffix         = fastq_out ? 'fastq' : 'fasta'
    def fastq_flag = fastq_out ? '--fastq-output' : ''
    def include_children = task.ext.include_children != null ? task.ext.include_children : true
    def include_parents  = task.ext.include_parents  != null ? task.ext.include_parents  : false
    def children_flag = include_children ? '--include-children' : ''
    def parents_flag  = include_parents  ? '--include-parents'  : ''
    def max_reads     = task.ext.max_reads ? "--max ${task.ext.max_reads}" : ''

    // ── Validate inputs ───────────────────────────────────────────────────
    if (!taxids || (taxids instanceof List && taxids.isEmpty())) {
        error("KRAKENTOOLS_EXTRACTKRAKENREADS: 'taxids' input is empty for sample '${meta.id}'")
    }
    def taxid_list = (taxids instanceof List) ? taxids.join(' ') : taxids.toString()

    // ── Resolve paired vs single end ──────────────────────────────────────
    def is_paired = reads instanceof List && reads.size() == 2
    def r1 = is_paired ? reads[0] : (reads instanceof List ? reads[0] : reads)
    def r2 = is_paired ? reads[1] : ''
    def s2_arg = is_paired ? "-s2 ${r2}" : ''
    def o2_arg = is_paired ? "-o2 ${prefix}.extracted_2.${suffix}" : ''

    """
    set -euo pipefail

    echo "[INFO] Sample          : ${meta.id}"
    echo "[INFO] Paired-end      : ${is_paired}"
    echo "[INFO] TaxIDs requested: ${taxid_list}"
    echo "[INFO] Include children: ${include_children}"
    echo "[INFO] Include parents : ${include_parents}"
    echo "[INFO] Output format   : ${suffix}"

    # ── Run extraction ───────────────────────────────────────────────────
    extract_kraken_reads.py \\
        -k ${kraken_output} \\
        -r ${kraken_report} \\
        -s ${r1} \\
        ${s2_arg} \\
        -o ${prefix}.extracted_1.${suffix} \\
        ${o2_arg} \\
        -t ${taxid_list} \\
        ${children_flag} \\
        ${parents_flag} \\
        ${fastq_flag} \\
        ${max_reads} \\
        ${args}

    # ── Sanity check: did we get any reads? ──────────────────────────────
    if [[ ! -s "${prefix}.extracted_1.${suffix}" ]]; then
        echo "[WARN] No reads extracted for taxids: ${taxid_list}" >&2
        # Create empty file with header so downstream doesn't crash
        touch "${prefix}.extracted_1.${suffix}"
    fi

    # ── Count extracted reads ────────────────────────────────────────────
    if [[ "${suffix}" == "fasta" ]]; then
        n_r1=\$(grep -c '^>' ${prefix}.extracted_1.${suffix} || echo 0)
    else
        n_r1=\$(( \$(wc -l < ${prefix}.extracted_1.${suffix}) / 4 ))
    fi

    n_r2=0
    if [[ "${is_paired}" == "true" && -s "${prefix}.extracted_2.${suffix}" ]]; then
        if [[ "${suffix}" == "fasta" ]]; then
            n_r2=\$(grep -c '^>' ${prefix}.extracted_2.${suffix} || echo 0)
        else
            n_r2=\$(( \$(wc -l < ${prefix}.extracted_2.${suffix}) / 4 ))
        fi
    fi

    # ── Write summary ────────────────────────────────────────────────────
    {
        echo "sample\\ttaxids\\tinclude_children\\tinclude_parents\\tn_extracted_R1\\tn_extracted_R2"
        echo "${meta.id}\\t${taxid_list}\\t${include_children}\\t${include_parents}\\t\${n_r1}\\t\${n_r2}"
    } > ${prefix}.summary.txt

    echo "[INFO] Extracted reads: R1=\${n_r1} R2=\${n_r2}"

    # ── Compress outputs ─────────────────────────────────────────────────
    if command -v pigz &>/dev/null; then
        pigz -p ${task.cpus} -f ${prefix}.extracted_1.${suffix}
    else
        gzip -f ${prefix}.extracted_1.${suffix}
    fi
    if [[ "${is_paired}" == "true" && -f "${prefix}.extracted_2.${suffix}" ]]; then
        if command -v pigz &>/dev/null; then
            pigz -p ${task.cpus} -f ${prefix}.extracted_2.${suffix}
        else
            gzip -f ${prefix}.extracted_2.${suffix}
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakentools: \$(extract_kraken_reads.py --version 2>&1 | sed 's/.*v//; s/ .*//' || echo "1.2")
        compression: \$(command -v pigz &>/dev/null && pigz --version 2>&1 | sed 's/pigz //' || gzip --version 2>&1 | head -1 | sed 's/gzip //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = (task.ext.fastq_output ?: false) ? 'fastq' : 'fasta'
    def is_paired = reads instanceof List && reads.size() == 2
    """
    echo ">stub_read1" > ${prefix}.extracted_1.${suffix}
    echo "ACGTACGTACGT" >> ${prefix}.extracted_1.${suffix}
    command -v pigz &>/dev/null && pigz -f ${prefix}.extracted_1.${suffix} || gzip -f ${prefix}.extracted_1.${suffix}
    # pigz -f ${prefix}.extracted_1.${suffix}

    if [[ "${is_paired}" == "true" ]]; then
        echo ">stub_read1" > ${prefix}.extracted_2.${suffix}
        echo "ACGTACGTACGT" >> ${prefix}.extracted_2.${suffix}
        command -v pigz &>/dev/null && pigz -f ${prefix}.extracted_2.${suffix} || gzip -f ${prefix}.extracted_2.${suffix}
        # pigz -f ${prefix}.extracted_2.${suffix}
    fi

    echo "sample\\ttaxids\\tn_extracted_R1" > ${prefix}.summary.txt
    echo "${meta.id}\\tstub\\t1" >> ${prefix}.summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakentools: 1.2
        compression: stub
    END_VERSIONS
    """
}

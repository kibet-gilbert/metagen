// =============================================================================
// modules/local/krakentools/extractkrakenreads/main.nf
//
// Extract reads assigned to target taxids from Kraken2 classification,
// with per-taxid random subsampling.
//
// Key behaviour:
//   - extract_kraken_reads.py is called ONCE PER TAXID (not once for all)
//     so that --max applies independently to each taxonomic group
//   - seqtk sample is applied per-taxid for true random subsampling when
//     the extracted count exceeds max_reads_per_taxid
//   - Per-taxid FASTAs are concatenated into final R1/R2 outputs
//   - A per-taxid extraction summary TSV is written alongside the global one
//
// Inputs:
//   tuple val(meta), path(reads), path(kraken_output), path(kraken_report)
//   val taxids   — list of NCBI taxids (integers or strings)
//
// Outputs:
//   extracted_r1     : R1 FASTA/FASTQ (all taxids combined)
//   extracted_r2     : R2 FASTA/FASTQ (all taxids combined, paired)
//   per_taxid_r1     : per-taxid R1 FASTAs (one file per taxid)
//   per_taxid_r2     : per-taxid R2 FASTAs (one file per taxid)
//   summary          : global extraction summary (compatible with VALIDATE_HITS)
//   per_taxid_summary: per-taxid counts before and after subsampling
//   versions         : versions.yml
//
// ext parameters (set in nextflow.config withName block):
//   ext.max_reads_per_taxid : int     — max reads per taxid after subsampling
//                                       (default: 50000)
//   ext.subsample_seed      : int     — seqtk random seed for reproducibility
//                                       (default: 42)
//   ext.include_children    : boolean — include descendant taxa (default: true)
//   ext.include_parents     : boolean — include ancestor taxa (default: false)
//   ext.fastq_output        : boolean — output FASTQ instead of FASTA (default: false)
//   ext.args                : string  — extra args passed to extract_kraken_reads.py
// =============================================================================

process KRAKENTOOLS_EXTRACTKRAKENREADS {
    tag "${meta.id}"
    label 'process_medium'

    // Both krakentools and seqtk needed — use a combined environment
    conda     "bioconda::krakentools=1.2 bioconda::seqtk=1.4 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/kibetgilbert/krakentools:1.2.1' :
        'docker.io/kibetgilbert/krakentools:1.2.1' }"
        // 'https://depot.galaxyproject.org/singularity/mulled-v2-0cdd3f94c78c03c0b99f4e93ba3810a5e62da4e3:13c4cc5c97f6ad4d0c2ceee6dd52e9f09b81b3b7-0' }"
        // 'biocontainers/mulled-v2-0cdd3f94c78c03c0b99f4e93ba3810a5e62da4e3:13c4cc5c97f6ad4d0c2ceee6dd52e9f09b81b3b7-0' }"
    // NOTE: the mulled container bundles krakentools + seqtk + pigz.
    // If this specific hash is unavailable on your HPC, build it locally:
    //   conda create -n krakentools_seqtk krakentools=1.2 seqtk=1.4 pigz=2.8
    //   apptainer build krakentools_seqtk.sif docker://quay.io/biocontainers/...
    // Or set ext.conda and ext.container per-process in your HPC config.

    publishDir(
        path:    { "${params.outdir}/extracted_reads/${meta.id}" },
        mode:    params.publish_dir_mode ?: 'copy',
        pattern: "*.{fasta.gz,fastq.gz,summary.txt,per_taxid_summary.tsv}",
        enabled: params.save_extracted_reads ?: true
    )

    input:
    tuple val(meta), path(reads), path(kraken_output), path(kraken_report)
    val   taxids

    output:
    tuple val(meta), path("${prefix}.extracted_1.${suffix}.gz"),                    emit: extracted_r1
    tuple val(meta), path("${prefix}.extracted_2.${suffix}.gz"),  optional: true,   emit: extracted_r2
    tuple val(meta), path("per_taxid/${prefix}_*.${suffix}.gz"),  optional: true,   emit: per_taxid_r1
    tuple val(meta), path("${prefix}.summary.txt"),                                  emit: summary
    tuple val(meta), path("${prefix}.per_taxid_summary.tsv"),                       emit: per_taxid_summary
    path  "versions.yml",                                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

script:
    prefix            = task.ext.prefix              ?: "${meta.id}"
    def args          = task.ext.args                ?: ''
    def fastq_out     = task.ext.fastq_output        ?: false
    suffix            = fastq_out ? 'fastq' : 'fasta'
    def fastq_flag    = fastq_out ? '--fastq-output' : ''
    def children_flag = (task.ext.include_children   ?: true)  ? '--include-children' : ''
    def parents_flag  = (task.ext.include_parents    ?: false) ? '--include-parents'  : ''
    def max_per_taxid = task.ext.max_reads_per_taxid ?: params.max_reads_per_taxid ?: 50000
    def seed          = task.ext.subsample_seed      ?: 42
    def is_paired     = reads instanceof List && reads.size() == 2
    def r1            = is_paired ? reads[0] : (reads instanceof List ? reads[0] : reads)
    def r2            = is_paired ? reads[1] : ''
    def s2_flag       = is_paired ? "-s2 ${r2}" : ''
    def taxid_str     = (taxids instanceof List)
                        ? taxids.collect { it.toString().trim() }.join(' ')
                        : taxids.toString().trim()
    """
    set -euo pipefail

    # ── Capture all Groovy compile-time values as named bash variables ─────
    # Groovy substitution happens ONCE here at the top.
    # Everything below uses \${BASH_VAR} (escaped) — evaluated at runtime.
    IS_PAIRED="${is_paired}"
    FORMAT="${suffix}"
    SAMPLE_PREFIX="${prefix}"
    MAX_READS=${max_per_taxid}
    SEED=${seed}
    CPUS=${task.cpus}

    echo "[INFO] Sample           : ${meta.id}"
    echo "[INFO] Paired-end       : \${IS_PAIRED}"
    echo "[INFO] Target taxids    : ${taxid_str}"
    echo "[INFO] Include children : ${children_flag}"
    echo "[INFO] Max per taxid    : \${MAX_READS}"
    echo "[INFO] Subsample seed   : \${SEED}"
    echo "[INFO] Output format    : \${FORMAT}"
    echo ""

    mkdir -p per_taxid
    : > \${SAMPLE_PREFIX}.extracted_1.\${FORMAT}
    if [[ "\${IS_PAIRED}" == "true" ]]; then
        : > \${SAMPLE_PREFIX}.extracted_2.\${FORMAT}
    fi

    printf "sample\ttaxid\tn_extracted_r1\tn_extracted_r2\tn_subsampled_r1\tn_subsampled_r2\tsubsampled\n" \
        > \${SAMPLE_PREFIX}.per_taxid_summary.tsv

    # ── Loop: one extraction call per taxid ───────────────────────────────
    for taxid in ${taxid_str}; do

        echo "[INFO] Processing taxid: \${taxid}"

        out_r1="per_taxid/\${SAMPLE_PREFIX}_\${taxid}_1.\${FORMAT}"
        out_r2="per_taxid/\${SAMPLE_PREFIX}_\${taxid}_2.\${FORMAT}"

        extract_kraken_reads.py \\
            -k ${kraken_output} \\
            -r ${kraken_report} \\
            -s  ${r1} \\
            ${s2_flag} \\
            -o  "\${out_r1}" \\
            ${is_paired ? "-o2 \"\${out_r2}\"" : ''} \\
            -t  \${taxid} \\
            ${children_flag} \\
            ${parents_flag} \\
            ${fastq_flag} \\
            ${args} \\
            2>&1 | grep -v '^\$' || true

        # Count extracted reads
        if [[ "\${FORMAT}" == "fasta" ]]; then
            n_r1=\$(grep -c '^>' "\${out_r1}" 2>/dev/null || echo 0)
        else
            n_r1=\$(awk 'NR%4==1' "\${out_r1}" 2>/dev/null | wc -l || echo 0)
        fi

        n_r2=0
        if [[ "\${IS_PAIRED}" == "true" && -s "\${out_r2}" ]]; then
            if [[ "\${FORMAT}" == "fasta" ]]; then
                n_r2=\$(grep -c '^>' "\${out_r2}" 2>/dev/null || echo 0)
            else
                n_r2=\$(awk 'NR%4==1' "\${out_r2}" 2>/dev/null | wc -l || echo 0)
            fi
        fi

        echo "[INFO]   taxid \${taxid}: extracted R1=\${n_r1} R2=\${n_r2}"

        if [[ \${n_r1} -eq 0 ]]; then
            echo "[INFO]   taxid \${taxid}: no reads extracted — skipping"
            printf "${meta.id}\t\${taxid}\t0\t0\t0\t0\tno_reads\n" \
                >> \${SAMPLE_PREFIX}.per_taxid_summary.tsv
            continue
        fi

        # ── seqtk random subsampling ─────────────────────────────────────
        # Correct argument order: seqtk sample -s SEED INPUT COUNT
        sub_r1="\${out_r1}"
        sub_r2="\${out_r2}"
        did_subsample="no"
        n_sub_r1=\${n_r1}
        n_sub_r2=\${n_r2}

        if [[ \${n_r1} -gt \${MAX_READS} ]]; then
            echo "[INFO]   taxid \${taxid}: subsampling \${n_r1} → \${MAX_READS} reads (seed=\${SEED})"
            sub_r1="per_taxid/\${SAMPLE_PREFIX}_\${taxid}_1.subsampled.\${FORMAT}"
            seqtk sample -s \${SEED} "\${out_r1}" \${MAX_READS} > "\${sub_r1}"
            n_sub_r1=\$(grep -c '^>' "\${sub_r1}" 2>/dev/null || \\
                        awk 'NR%4==1' "\${sub_r1}" 2>/dev/null | wc -l || echo 0)

            if [[ "\${IS_PAIRED}" == "true" && -s "\${out_r2}" ]]; then
                sub_r2="per_taxid/\${SAMPLE_PREFIX}_\${taxid}_2.subsampled.\${FORMAT}"
                seqtk sample -s \${SEED} "\${out_r2}" \${MAX_READS} > "\${sub_r2}"
                n_sub_r2=\$(grep -c '^>' "\${sub_r2}" 2>/dev/null || \\
                            awk 'NR%4==1' "\${sub_r2}" 2>/dev/null | wc -l || echo 0)
            fi
            did_subsample="yes"
        fi

        echo "[INFO]   taxid \${taxid}: final R1=\${n_sub_r1} R2=\${n_sub_r2} (subsampled=\${did_subsample})"

        printf "${meta.id}\t\${taxid}\t\${n_r1}\t\${n_r2}\t\${n_sub_r1}\t\${n_sub_r2}\t\${did_subsample}\n" \
            >> \${SAMPLE_PREFIX}.per_taxid_summary.tsv

        cat "\${sub_r1}" >> \${SAMPLE_PREFIX}.extracted_1.\${FORMAT}
        if [[ "\${IS_PAIRED}" == "true" && -s "\${sub_r2}" ]]; then
            cat "\${sub_r2}" >> \${SAMPLE_PREFIX}.extracted_2.\${FORMAT}
        fi

    done

    # ── Totals ───────────────────────────────────────────────────────────
    if [[ "\${FORMAT}" == "fasta" ]]; then
        total_r1=\$(grep -c '^>' \${SAMPLE_PREFIX}.extracted_1.\${FORMAT} 2>/dev/null || echo 0)
    else
        total_r1=\$(awk 'NR%4==1' \${SAMPLE_PREFIX}.extracted_1.\${FORMAT} 2>/dev/null | wc -l || echo 0)
    fi
    total_r2=0
    if [[ "\${IS_PAIRED}" == "true" && -s "\${SAMPLE_PREFIX}.extracted_2.\${FORMAT}" ]]; then
        if [[ "\${FORMAT}" == "fasta" ]]; then
            total_r2=\$(grep -c '^>' \${SAMPLE_PREFIX}.extracted_2.\${FORMAT} 2>/dev/null || echo 0)
        else
            total_r2=\$(awk 'NR%4==1' \${SAMPLE_PREFIX}.extracted_2.\${FORMAT} 2>/dev/null | wc -l || echo 0)
        fi
    fi

    echo "[INFO] Total combined R1=\${total_r1} R2=\${total_r2}"

    {
        printf "sample\ttaxids\tinclude_children\tmax_per_taxid\tseed\tn_extracted_R1\tn_extracted_R2\n"
        printf "${meta.id}\t${taxid_str}\t${children_flag ?: 'no'}\t\${MAX_READS}\t\${SEED}\t\${total_r1}\t\${total_r2}\n"
    } > \${SAMPLE_PREFIX}.summary.txt

    echo "[INFO] Compressing outputs..."
    pigz -p \${CPUS} -f \${SAMPLE_PREFIX}.extracted_1.\${FORMAT}

    if [[ "\${IS_PAIRED}" == "true" && -f "\${SAMPLE_PREFIX}.extracted_2.\${FORMAT}" ]]; then
        pigz -p \${CPUS} -f \${SAMPLE_PREFIX}.extracted_2.\${FORMAT}
    fi

    for f in per_taxid/\${SAMPLE_PREFIX}_*.\${FORMAT}; do
        [[ -f "\${f}" ]] && pigz -p 2 -f "\${f}" || true
    done

    echo "[INFO] Done: ${meta.id}"
    cat \${SAMPLE_PREFIX}.per_taxid_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakentools: \$(extract_kraken_reads.py --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+' | head -1 || echo "1.2")
        seqtk: \$(seqtk 2>&1 | grep -oP 'Version: \\K[^ ]+' || echo "1.4")
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //')
    END_VERSIONS
    """
}

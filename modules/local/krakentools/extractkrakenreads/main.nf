// =============================================================================
// modules/local/krakentools/extractkrakenreads/main.nf
//
// Extract and subsample reads by taxid from Kraken2 results.
// Python logic lives in bin/extract_reads_by_taxid.py (auto-added to PATH).
//
// Sub-steps:
//   1a. parse-report    -> taxonomy tree, abundance table, taxid groups
//   1b. extract-readids -> read_id/taxid pairs from kraken2 output
//   1c. subsample       -> per-group random subsampling without replacement
//   1d. seqtk subseq   -> physical read extraction from R1 / R2 FASTQ
// =============================================================================

process KRAKENTOOLS_EXTRACTKRAKENREADS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::seqtk=1.4 conda-forge::python=3.10 conda-forge::pigz"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/kibetgilbert/krakentools:1.2.1' :
        'docker.io/kibetgilbert/krakentools:1.2.1' }"

    publishDir(
        path:    { "${params.outdir}/extracted_reads/${meta.id}" },
        mode:    params.publish_dir_mode ?: 'copy',
        pattern: "*.{fasta.gz,fastq.gz,txt,tsv,json}",
        enabled: params.save_extracted_reads ?: true
    )

    input:
    tuple val(meta), path(reads), path(kraken_output), path(kraken_report)
    val   taxids

    output:
    tuple val(meta), path("${prefix}.extracted_R1.fasta.gz"),   emit: extracted_r1
    tuple val(meta), path("${prefix}.extracted_R2.fasta.gz"),   emit: extracted_r2,          optional: true
    tuple val(meta), path("${prefix}.taxid_abundance.tsv"),     emit: taxid_abundance
    tuple val(meta), path("${prefix}.subsampling_summary.tsv"), emit: subsampling_summary
    tuple val(meta), path("${prefix}.readid_taxid.tsv"),        emit: readid_taxid
    tuple val(meta), path("${prefix}.taxid_groups.json"),       emit: taxid_groups
    tuple val(meta), path("${prefix}.summary.txt"),             emit: summary
    path  "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix        = task.ext.prefix              ?: "${meta.id}"
    def max_reads = task.ext.max_reads_per_taxid ?: params.max_reads_per_taxid ?: 50000
    def seed      = task.ext.subsample_seed      ?: 42
    def is_paired = reads instanceof List && reads.size() == 2
    def r1        = is_paired ? reads[0] : reads
    def r2        = is_paired ? reads[1] : 'NO_R2'
    def taxid_csv = (taxids instanceof List)
                    ? taxids.collect { it.toString().trim() }.join(',')
                    : taxids.toString().trim()
    """
    set -euo pipefail

    # Resolve Groovy values into named bash variables ONCE.
    # Everything below uses \${BASH_VAR} (escaped dollar, evaluated at runtime).
    IS_PAIRED="${is_paired}"
    R1="${r1}"
    R2="${r2}"
    PREFIX="${prefix}"
    MAX_READS="${max_reads}"
    SEED="${seed}"
    TAXID_CSV="${taxid_csv}"
    CPUS="${task.cpus}"

    echo "========================================================"
    echo "[INFO] Sample       : ${meta.id}"
    echo "[INFO] Taxids       : \${TAXID_CSV}"
    echo "[INFO] Max reads    : \${MAX_READS} per taxid group"
    echo "[INFO] Seed         : \${SEED}"
    echo "[INFO] Paired       : \${IS_PAIRED}"
    echo "========================================================"

    mkdir -p per_taxid

    # ── STEP 1a: Parse kraken2 report ──────────────────────────────────────
    echo "[STEP 1a] Parsing kraken2 report..."
    extract_reads_by_taxid.py parse-report \\
        --report  ${kraken_report} \\
        --taxids  "\${TAXID_CSV}" \\
        --prefix  "\${PREFIX}"

    echo "[STEP 1a] Done. Taxids found: \$(wc -l < \${PREFIX}.all_taxids.txt)"

    # ── STEP 1b: Extract read_id/taxid pairs ───────────────────────────────
    echo "[STEP 1b] Extracting read IDs from kraken2 output..."
    extract_reads_by_taxid.py extract-readids \\
        --output  ${kraken_output} \\
        --taxids  "\${PREFIX}.all_taxids.txt" \\
        --prefix  "\${PREFIX}"

    n_readids=\$(( \$(wc -l < \${PREFIX}.readid_taxid.tsv) - 1 ))
    echo "[STEP 1b] Done. Read IDs extracted: \${n_readids}"

    # ── STEP 1c: Subsample per taxid group ─────────────────────────────────
    echo "[STEP 1c] Subsampling read IDs per taxid group..."
    extract_reads_by_taxid.py subsample \\
        --readids   "\${PREFIX}.readid_taxid.tsv" \\
        --groups    "\${PREFIX}.taxid_groups.json" \\
        --max-reads "\${MAX_READS}" \\
        --seed      "\${SEED}" \\
        --prefix    "\${PREFIX}"

    echo "[STEP 1c] Done."
    cat \${PREFIX}.subsampling_summary.tsv

    # ── STEP 1d: Extract reads from FASTQ by read ID list ──────────────────
    echo "[STEP 1d] Extracting reads with seqtk subseq..."

    : > \${PREFIX}.extracted_R1.fasta
    if [[ "\${IS_PAIRED}" == "true" ]]; then
        : > \${PREFIX}.extracted_R2.fasta
    fi

    # Auto-detect mate suffix style from the actual FASTQ headers.
    # Handles /1 /2, .1 .2, or no suffix at all — avoids hardcoding.
    R1_HEADER=\$(gzip -cd "\${R1}" | sed -n '1p')
    R1_SUFFIX=""
    if   [[ "\${R1_HEADER}" == *"/1" ]]; then R1_SUFFIX="/1"
    elif [[ "\${R1_HEADER}" == *".1" ]]; then R1_SUFFIX=".1"
    fi
    echo "[INFO] Detected R1 suffix: '\${R1_SUFFIX}'"

    R2_SUFFIX=""
    if [[ "\${IS_PAIRED}" == "true" ]]; then
        R2_HEADER=\$(gzip -cd "\${R2}" | sed -n '1p')
        if   [[ "\${R2_HEADER}" == *"/2" ]]; then R2_SUFFIX="/2"
        elif [[ "\${R2_HEADER}" == *".2" ]]; then R2_SUFFIX=".2"
        fi
        echo "[INFO] Detected R2 suffix: '\${R2_SUFFIX}'"
    fi

    for id_file in per_taxid/\${PREFIX}_*.readids.txt; do
        [[ -f "\${id_file}" ]] || continue

        n_ids=\$(wc -l < "\${id_file}" || echo 0)
        taxid=\$(basename "\${id_file}" .readids.txt | sed "s/\${PREFIX}_//")

        if [[ \${n_ids} -eq 0 ]]; then
            echo "[INFO]   taxid \${taxid}: 0 read IDs — skipping"
            continue
        fi

        echo "[INFO]   taxid \${taxid}: extracting \${n_ids} read pairs"

        # Build R1-specific and R2-specific ID lists with correct suffix
        r1_ids="per_taxid/\${PREFIX}_\${taxid}.r1ids.txt"
        sed "s|\$|\${R1_SUFFIX}|" "\${id_file}" > "\${r1_ids}"
        seqtk subseq "\${R1}" "\${r1_ids}" >> \${PREFIX}.extracted_R1.fasta

        if [[ "\${IS_PAIRED}" == "true" ]]; then
            r2_ids="per_taxid/\${PREFIX}_\${taxid}.r2ids.txt"
            sed "s|\$|\${R2_SUFFIX}|" "\${id_file}" > "\${r2_ids}"
            seqtk subseq "\${R2}" "\${r2_ids}" >> \${PREFIX}.extracted_R2.fasta
        fi

        n_r1=\$(awk '/^>/{c++} END{print c+0}' \${PREFIX}.extracted_R1.fasta)
        # n_r1=\$(grep -c '^>' \${PREFIX}.extracted_R1.fasta 2>/dev/null || echo 0)
        echo "[INFO]   taxid \${taxid}: R1 sequences so far: \${n_r1}"
    done

    # ── Final counts ────────────────────────────────────────────────────────
    total_r1=\$(grep -c '^>' \${PREFIX}.extracted_R1.fasta 2>/dev/null || echo 0)
    total_r2=0
    if [[ "\${IS_PAIRED}" == "true" && -s "\${PREFIX}.extracted_R2.fasta" ]]; then
        total_r2=\$(grep -c '^>' \${PREFIX}.extracted_R2.fasta 2>/dev/null || echo 0)
    fi

    echo ""
    echo "[INFO] Total extracted: R1=\${total_r1}  R2=\${total_r2}"

    # ── Global summary for VALIDATE_HITS ────────────────────────────────────
    {
        printf "sample\ttaxids\tmax_per_taxid\tseed\tn_extracted_R1\tn_extracted_R2\n"
        printf "${meta.id}\t\${TAXID_CSV}\t\${MAX_READS}\t\${SEED}\t\${total_r1}\t\${total_r2}\n"
    } > \${PREFIX}.summary.txt

    # ── Compress ─────────────────────────────────────────────────────────────
    pigz -p \${CPUS} -f \${PREFIX}.extracted_R1.fasta
    if [[ "\${IS_PAIRED}" == "true" && -f "\${PREFIX}.extracted_R2.fasta" ]]; then
        pigz -p \${CPUS} -f \${PREFIX}.extracted_R2.fasta
    fi

    echo "[INFO] Done: ${meta.id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        seqtk: \$(seqtk 2>&1 | grep -oP 'Version: \\K[0-9.]+' || echo "1.4")
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p per_taxid
    printf ">stub\\nACGT\\n" | pigz > ${prefix}.extracted_R1.fasta.gz
    printf ">stub\\nACGT\\n" | pigz > ${prefix}.extracted_R2.fasta.gz
    printf "target_taxid\ttaxid\n"           > ${prefix}.taxid_abundance.tsv
    printf "target_taxid\tn_reads_sampled\n" > ${prefix}.subsampling_summary.tsv
    printf "read_id\ttaxid\n"               > ${prefix}.readid_taxid.tsv
    printf "{}\n"                           > ${prefix}.taxid_groups.json
    printf "sample\tn_extracted_R1\n"       > ${prefix}.summary.txt

    {
    echo '"${task.process}":'
    echo "    python: \$(python3 --version | sed 's/Python //')"
    echo "    seqtk: \$(seqtk 2>&1 | grep -oP 'Version: \\K[0-9.]+' || echo '1.4')"
    echo "    pigz: \$(pigz --version 2>&1 | sed 's/pigz //')"
    } > versions.yml
    """
}

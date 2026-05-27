// =============================================================================
// modules/local/blast/blastn/main.nf
//
// BLASTn alignment for pathogen validation.
//
// Supports three database modes via task.ext.db_mode:
//   1. 'local'   : local BLAST DB on shared filesystem (e.g. NCBI nt downloaded to HPC)
//   2. 'remote'  : NCBI remote BLAST (-remote flag) — slow, rate-limited, internet required
//   3. 'custom'  : custom local BLAST DB (e.g. curated pathogen genomes)
//
// Input format is FASTA (one or more sequences per file).
// Designed to accept output from KRAKENTOOLS_EXTRACTKRAKENREADS.
//
// Inputs:
//   tuple val(meta), path(fasta)  : FASTA file(s) to BLAST
//   path  db_dir                   : BLAST DB directory (when db_mode='local'/'custom'; pass '[]' for remote)
//   val   db_name                  : DB name (e.g. 'nt', 'core_nt', 'refseq_genomic')
//
// Outputs:
//   tuple val(meta), path("*.blastn.tsv")     : BLAST tabular results (outfmt 6/7)
//   tuple val(meta), path("*.blastn.summary.tsv") : per-taxid hit summary
//   path  "versions.yml"
// =============================================================================

process BLAST_BLASTN {
    tag "${meta.id}_${db_name}"
    label 'process_high'

    conda     "bioconda::blast=2.15.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1' :
        'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    publishDir(
        path:  { "${params.outdir}/blast/${db_name}/${meta.id}" },
        mode:  params.publish_dir_mode ?: 'copy',
        pattern: "*.{tsv,log}"
    )

    input:
    tuple val(meta), path(fasta)
    path  db_dir
    val   db_name

    output:
    tuple val(meta), path("${prefix}.blastn.tsv"),                    emit: tsv
    tuple val(meta), path("${prefix}.blastn.summary.tsv"),            emit: summary
    tuple val(meta), path("${prefix}.blastn.log"),                    emit: log
    path  "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ── Resolve parameters ────────────────────────────────────────────────
    prefix             = task.ext.prefix          ?: "${meta.id}__${db_name}"
    def args           = task.ext.args            ?: ''
    def db_mode        = task.ext.db_mode         ?: 'local'    // local | remote | custom
    def evalue         = task.ext.evalue          ?: '1e-10'
    def max_target_seqs= task.ext.max_target_seqs ?: 5
    def perc_identity  = task.ext.perc_identity   ?: 90
    def qcov_hsp       = task.ext.qcov_hsp        ?: 80
    def num_alignments = task.ext.num_alignments  ?: 5
    def word_size      = task.ext.word_size       ?: 28        // megablast default; lower for divergent
    def task_mode      = task.ext.task_mode       ?: 'megablast' // megablast|dc-megablast|blastn|blastn-short

    // ── Output format: tabular with extended fields including taxonomy ───
    def outfmt = task.ext.outfmt ?: '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames qcovs qcovhsp stitle'
    def outfmt_header = (outfmt =~ /^[6-7] (.*)/)[0][1]?.replaceAll(' ', '\\t') ?: ''

    // ── Mode-specific flags ───────────────────────────────────────────────
    def db_flag
    def remote_flag
    def threads_flag
    if (db_mode == 'remote') {
        db_flag      = "-db ${db_name}"
        remote_flag  = '-remote'
        threads_flag = ''   // -num_threads cannot be used with -remote
    } else {
        db_flag      = "-db ${db_dir}/${db_name}"
        remote_flag  = ''
        threads_flag = "-num_threads ${task.cpus}"
    }

    """
    set -euo pipefail

    echo "[INFO] Sample          : ${meta.id}"
    echo "[INFO] Database        : ${db_name} (${db_mode})"
    echo "[INFO] Task mode       : ${task_mode}"
    echo "[INFO] E-value         : ${evalue}"
    echo "[INFO] Min identity    : ${perc_identity}%"
    echo "[INFO] Min query cov   : ${qcov_hsp}%"

    # ── Verify input has sequences ────────────────────────────────────────
    if [[ "${fasta}" == *.gz ]]; then
        n_seqs=\$(zgrep -c '^>' ${fasta} || echo 0)
        # Decompress for blast (some versions don't read gz)
        zcat ${fasta} > input.fa
        FASTA_INPUT="input.fa"
    else
        n_seqs=\$(grep -c '^>' ${fasta} || echo 0)
        FASTA_INPUT="${fasta}"
    fi

    echo "[INFO] Input sequences : \${n_seqs}"

    if [[ "\${n_seqs}" -eq 0 ]]; then
        echo "[WARN] No input sequences — emitting empty result files" >&2
        # Write header row so downstream parsers don't break
        echo -e "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tstaxids\\tsscinames\\tscomnames\\tqcovs\\tqcovhsp\\tstitle" > ${prefix}.blastn.tsv
        echo -e "taxid\\tsci_name\\tn_hits\\tmean_pident\\tmean_qcovs\\tmean_bitscore\\tbest_evalue" > ${prefix}.blastn.summary.tsv
        echo "no_input" > ${prefix}.blastn.log

        cat <<EOF > versions.yml
    "${task.process}":
        blast: \$(blastn -version | head -1 | sed 's/blastn: //')
    EOF
        exit 0
    fi

    # ── Run BLASTn ────────────────────────────────────────────────────────
    blastn \\
        -task ${task_mode} \\
        -query \${FASTA_INPUT} \\
        ${db_flag} \\
        ${remote_flag} \\
        ${threads_flag} \\
        -evalue ${evalue} \\
        -perc_identity ${perc_identity} \\
        -qcov_hsp_perc ${qcov_hsp} \\
        -word_size ${word_size} \\
        -max_target_seqs ${max_target_seqs} \\
        -outfmt "${outfmt}" \\
        ${args} \\
        > ${prefix}.blastn.tsv \\
        2> ${prefix}.blastn.log

    # ── Prepend header to TSV (BLAST outfmt 6 has no header by default) ──
    if [[ -s "${prefix}.blastn.tsv" ]]; then
        # Build header from outfmt definition
        header=\$(echo "${outfmt}" | sed -E 's/^[6-7] //; s/ /\\t/g')
        # sed -i "1i\\\${header}" ${prefix}.blastn.tsv
        { echo "\${header}"; cat ${prefix}.blastn.tsv; } > ${prefix}.blastn.tmp
        mv ${prefix}.blastn.tmp ${prefix}.blastn.tsv
    fi

    n_hits=\$(( \$(wc -l < ${prefix}.blastn.tsv) - 1 ))
    echo "[INFO] BLAST hits      : \${n_hits}"

    # ── Build per-taxid summary ──────────────────────────────────────────
    awk -F'\\t' '
    NR==1 {
        for(i=1;i<=NF;i++) col[\$i] = i
        next
    }
    NR>1 && \$col["staxids"] != "" {
        # Handle multi-taxid hits (semicolon-separated)
        n = split(\$col["staxids"], taxa, ";")
        for(i=1;i<=n;i++) {
            t = taxa[i]
            hits[t]++
            sum_pident[t]   += \$col["pident"]
            sum_qcovs[t]    += \$col["qcovs"]
            sum_bitscore[t] += \$col["bitscore"]
            if(best_eval[t]=="" || \$col["evalue"]+0 < best_eval[t]+0) best_eval[t] = \$col["evalue"]
            sci_name[t] = \$col["sscinames"]
        }
    }
    END {
        print "taxid\\tsci_name\\tn_hits\\tmean_pident\\tmean_qcovs\\tmean_bitscore\\tbest_evalue"
        for(t in hits) {
            printf "%s\\t%s\\t%d\\t%.2f\\t%.2f\\t%.2f\\t%s\\n",
                t, sci_name[t], hits[t],
                sum_pident[t]/hits[t],
                sum_qcovs[t]/hits[t],
                sum_bitscore[t]/hits[t],
                best_eval[t]
        }
    }' ${prefix}.blastn.tsv \\
    | (head -1; tail -n +2 | sort -t\$'\\t' -k3,3nr) \\
    > ${prefix}.blastn.summary.tsv

    n_taxa=\$(( \$(wc -l < ${prefix}.blastn.summary.tsv) - 1 ))
    echo "[INFO] Unique taxa hit : \${n_taxa}"

    # Cleanup decompressed input
    [[ -f input.fa ]] && rm -f input.fa

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version | head -1 | sed 's/blastn: //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}__${db_name}"
    """
    echo -e "qseqid\\tsseqid\\tpident\\tlength\\tstaxids\\tsscinames" > ${prefix}.blastn.tsv
    echo -e "read1\\tref1\\t99.5\\t150\\t562\\tEscherichia coli" >> ${prefix}.blastn.tsv

    echo -e "taxid\\tsci_name\\tn_hits\\tmean_pident\\tmean_qcovs\\tmean_bitscore\\tbest_evalue" > ${prefix}.blastn.summary.tsv
    echo -e "562\\tEscherichia coli\\t1\\t99.50\\t100.00\\t280.00\\t1e-50" >> ${prefix}.blastn.summary.tsv

    touch ${prefix}.blastn.log

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        blast: 2.15.0
    END_VERSIONS
    """
}

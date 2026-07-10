// =============================================================================
// PROCESS: Cross-validate Kraken2 taxid claims vs BLAST hit evidence
// =============================================================================
process VALIDATE_HITS {
    tag "${meta.id}"
    label 'process_low'

    conda     "conda-forge::python=3.14 conda-forge::pandas=2.2.1"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    publishDir(
        path:    "${params.outdir}/validation",
        mode:    params.publish_dir_mode ?: 'copy',
        pattern: "*.{tsv,json}"
    )

    input:
    tuple val(meta), path(kraken_summary), path(blast_summary)
    val   target_taxids

    output:
    tuple val(meta), path("${meta.id}.validation.tsv"),  emit: report
    tuple val(meta), path("${meta.id}.validation.json"), emit: json
    path  "versions.yml",                                  emit: versions

    script:
    def taxid_str    = (target_taxids instanceof List) ? target_taxids.join(',') : target_taxids.toString()
    def min_hits     = task.ext.min_blast_hits ?: params.min_blast_hits ?: 3
    def min_pident   = task.ext.min_pident     ?: params.min_pident     ?: 95.0
    def min_qcovs    = task.ext.min_qcovs      ?: params.min_qcovs      ?: 80.0
    def min_bitscore = task.ext.min_bitscore   ?: params.min_bitscore   ?: 100.0
    def has_blast    = blast_summary instanceof List ? 'False' : 'True'
    """
    #!/usr/bin/env python3
    import pandas as pd, json, sys, os

    SAMPLE       = "${meta.id}"
    TARGET_STR   = "${taxid_str}"
    MIN_HITS     = ${min_hits}
    MIN_PIDENT   = ${min_pident}
    MIN_QCOVS    = ${min_qcovs}
    MIN_BITSCORE = ${min_bitscore}
    HAS_BLAST    = ${has_blast}

    targets = [t.strip() for t in TARGET_STR.split(',') if t.strip()]

    blast_path = "${blast_summary}"
    if HAS_BLAST and os.path.isfile(blast_path) and os.path.getsize(blast_path) > 0:
        try:
            blast_df = pd.read_csv(blast_path, sep='\\t', dtype={'taxid': str})
        except Exception:
            blast_df = pd.DataFrame()
    else:
        blast_df = pd.DataFrame()

    try:
        kraken_df    = pd.read_csv("${kraken_summary}", sep='\\t')
        reads_r1     = int(kraken_df.get('n_extracted_R1', [0]).iloc[0])
        reads_r2     = int(kraken_df.get('n_extracted_R2', [0]).iloc[0])
    except Exception:
        reads_r1 = reads_r2 = 0

    records = []
    for tax in targets:
        base = {'sample': SAMPLE, 'target_taxid': tax,
                'kraken_reads_r1': reads_r1, 'kraken_reads_r2': reads_r2}

        if blast_df.empty or 'taxid' not in blast_df.columns:
            records.append({**base, 'sci_name': 'NA', 'blast_hits': 0,
                'blast_mean_pident': None, 'blast_mean_qcovs': None,
                'blast_mean_bitscore': None, 'best_evalue': None,
                'validation_status': 'NOT_VALIDATED', 'checks_passed': '0/4',
                'reason': 'No BLAST output'})
            continue

        hit = blast_df[blast_df['taxid'] == tax]
        if hit.empty:
            n_other = len(blast_df)
            records.append({**base, 'sci_name': 'NA', 'blast_hits': 0,
                'blast_mean_pident': None, 'blast_mean_qcovs': None,
                'blast_mean_bitscore': None, 'best_evalue': None,
                'validation_status': 'FAILED', 'checks_passed': '0/4',
                'reason': f'No BLAST hits to taxid {tax}; {n_other} hits to other taxa'})
        else:
            n  = int(hit['n_hits'].iloc[0])
            pi = float(hit['mean_pident'].iloc[0])
            qc = float(hit['mean_qcovs'].iloc[0])
            bs = float(hit['mean_bitscore'].iloc[0])
            ev = hit['best_evalue'].iloc[0]
            sc = hit['sci_name'].iloc[0]
            passed = sum([n >= MIN_HITS, pi >= MIN_PIDENT, qc >= MIN_QCOVS, bs >= MIN_BITSCORE])
            if   passed == 4: status, reason = 'CONFIRMED',     'All 4 thresholds passed'
            elif passed >= 2: status, reason = 'WEAK_EVIDENCE', f'{passed}/4 thresholds passed'
            else:             status, reason = 'FAILED',        f'Only {passed}/4 thresholds passed'
            records.append({**base, 'sci_name': sc, 'blast_hits': n,
                'blast_mean_pident': round(pi,2), 'blast_mean_qcovs': round(qc,2),
                'blast_mean_bitscore': round(bs,2), 'best_evalue': ev,
                'validation_status': status, 'checks_passed': f'{passed}/4', 'reason': reason})

    pd.DataFrame(records).to_csv(f'{SAMPLE}.validation.tsv', sep='\\t', index=False)
    with open(f'{SAMPLE}.validation.json', 'w') as fh:
        json.dump(records, fh, indent=2, default=str)
    print(f'[INFO] {SAMPLE}: wrote {len(records)} records')

    with open('versions.yml', 'w') as fh:
        fh.write('"${task.process}":\\n')
        fh.write(f'    python: {sys.version.split()[0]}\\n')
        fh.write(f'    pandas: {pd.__version__}\\n')
    """

    stub:
    """
    printf "sample\\ttarget_taxid\\tvalidation_status\\n${meta.id}\\tstub\\tCONFIRMED\\n" \
      > ${meta.id}.validation.tsv
    printf '[]\\n' > ${meta.id}.validation.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10
        pandas: 2.0
    END_VERSIONS
    """
}



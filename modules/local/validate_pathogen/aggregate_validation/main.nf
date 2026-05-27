// =============================================================================
// PROCESS: Aggregate per-sample validation TSVs into one master report
// =============================================================================
process AGGREGATE_VALIDATION {
    label 'process_low'

    conda     "conda-forge::python=3.10 conda-forge::pandas=2.0"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.3' :
        'quay.io/biocontainers/pandas:2.0.3' }"

    publishDir(
        path:    "${params.outdir}/validation_summary",
        mode:    params.publish_dir_mode ?: 'copy',
        pattern: "*.{tsv,txt}"
    )

    input:
    path(tsvs)   // collected list of *.validation.tsv files

    output:
    path "all_samples_validation.tsv", emit: tsv
    path "validation_summary.txt",     emit: summary
    path "versions.yml",               emit: versions

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd, glob, sys

    files = sorted(glob.glob('*.validation.tsv'))
    print(f'[INFO] Aggregating {len(files)} files')

    dfs = []
    for f in files:
        try:
            dfs.append(pd.read_csv(f, sep='\\t'))
        except Exception as e:
            print(f'[WARN] Cannot read {f}: {e}', file=sys.stderr)

    if not dfs:
        open('all_samples_validation.tsv', 'w').close()
        open('validation_summary.txt', 'w').close()
    else:
        combined = pd.concat(dfs, ignore_index=True)
        combined.sort_values(['target_taxid','sample']).to_csv(
            'all_samples_validation.tsv', sep='\\t', index=False)

        status_counts = combined.groupby('validation_status').size()
        per_taxon     = combined.groupby(['target_taxid','sci_name','validation_status']).size().reset_index(name='n')

        with open('validation_summary.txt','w') as fh:
            fh.write("=== Pathogen Validation Summary ===\\n\\n")
            fh.write(f"Sample x taxid pairs  : {len(combined)}\\n")
            fh.write(f"Unique samples        : {combined['sample'].nunique()}\\n")
            fh.write(f"Target taxids         : {combined['target_taxid'].nunique()}\\n\\n")
            fh.write("--- Status counts ---\\n")
            fh.write(status_counts.to_string())
            fh.write("\\n\\n--- Per-taxon outcomes ---\\n")
            fh.write(per_taxon.to_string(index=False))
            fh.write("\\n")
        print(status_counts.to_string())

    with open('versions.yml','w') as fh:
        fh.write('"${task.process}":\\n')
        fh.write(f'    python: {sys.version.split()[0]}\\n')
        fh.write(f'    pandas: {pd.__version__}\\n')
    """

    stub:
    """
    printf "sample\\ttarget_taxid\\tvalidation_status\\n" > all_samples_validation.tsv
    printf "=== Stub ===\\n"                               > validation_summary.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10
        pandas: 2.0
    END_VERSIONS
    """
}



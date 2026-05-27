/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/ccmetagen_merge.nf
    Merge per-sample CCMetagen CSV tables into a single combined abundance table.

    Uses CCMetagen's built-in merge utility (CCMetagen_merge.py) when available,
    with a fallback to a plain awk/Python merge.

    Input:
        tuple val(meta), path(csvs, stageAs: "?/*")  // collected per-sample CSVs
    Output:
        tuple val(meta), path("*.combined.csv")    , emit: combined
        path  "versions.yml"                       , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CCMETAGEN_MERGE {

    tag "all_samples"

    label 'process_single'

    // conda "bioconda::ccmetagen=1.1.4 conda-forge::pandas=2.2"
    // container "${ workflow.containerEngine == 'singularity' ?
    //     'oras://community.wave.seqera.io/library/ccmetagen:1.1.4' :
    //     'biocontainers/ccmetagen:1.1.4' }"

    input:
    tuple val(meta), path(csvs, stageAs: "?/*")

    output:
    tuple val(meta), path("*.combined.csv") , emit: combined
    path  "versions.yml"                    , emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: (meta.id ?: 'all_samples')

    """
    # ── Collect all CSV paths into a flat list ────────────────────────────
    find . -name "*.csv" -not -name "*.combined.csv" | sort > csv_list.txt
    echo "[CCMETAGEN_MERGE] Found \$(wc -l < csv_list.txt) sample CSVs"

    # ── Try CCMetagen_merge.py first (preferred) ──────────────────────────
    if command -v CCMetagen_merge.py &>/dev/null; then
        echo "Using CCMetagen_merge.py"
        CCMetagen_merge.py \\
            --input_fp . \\
            --output_fp ${prefix}.combined.csv \\
            ${args}

    else
        # ── Fallback: Python pandas merge ────────────────────────────────
        echo "CCMetagen_merge.py not found — using Python pandas fallback"
        python3 - <<'PYEOF'
import sys, pathlib, pandas as pd

csv_paths = sorted(pathlib.Path('.').rglob('*.csv'))
csv_paths = [p for p in csv_paths if '.combined' not in p.name]
print(f"Merging {len(csv_paths)} CSVs...")

frames = []
for p in csv_paths:
    try:
        df = pd.read_csv(p)
        sample = p.parent.name if p.parent.name != '.' else p.stem
        df.insert(0, 'sample', sample)
        frames.append(df)
    except Exception as e:
        print(f"WARNING: could not read {p}: {e}", file=sys.stderr)

if not frames:
    sys.exit("ERROR: no valid CSV files found to merge")

merged = pd.concat(frames, ignore_index=True)
out = "${prefix}.combined.csv"
merged.to_csv(out, index=False)
print(f"Written: {out}  ({len(merged)} rows, {merged.shape[1]} columns)")
PYEOF
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CCMetagen_merge: \$(CCMetagen_merge.py --version 2>&1 | head -1 || echo "fallback-pandas")
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: (meta.id ?: 'all_samples')
    """
    echo "sample,kingdom,phylum,class,order,family,genus,species,reads" \\
        > ${prefix}.combined.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CCMetagen_merge: 1.1.4
        python: 3.11.0
    END_VERSIONS
    """
}

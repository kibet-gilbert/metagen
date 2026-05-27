/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    modules/local/fasta_clean.nf
    MODULE 2 of 2 in the CCMetagen DB-build chain.

    Cleans a raw FASTA (from BLAST_EXPORT or download_all_refseq.sh) and builds
    a KMA index ready for CCMetagen alignment.

    Steps executed by clean_fasta.sh:
      1. Title filtering  — drops PREDICTED, vector, synthetic, uncultured …
      2. Length filter    — keeps sequences >= MIN_LEN (default 300)
      3. Deduplication    — seqkit rmdup by sequence  (optional, env DEDUP_SEQUENCES)
      4. Low-complexity   — dustmasker → hard-mask to N  (optional, env MASK_LOW_COMPLEX)
      5. Compress         — pigz
      6. KMA index        — kma index → output DB prefix

    Input:
        tuple val(meta), path(raw_fasta)     // .fasta or .fasta.gz from BLAST_EXPORT
    Output:
        tuple val(meta), path("*.ccmetagen.fasta.gz")      , emit: cleaned_fasta
        tuple val(meta), path("*_ccmetagen.{name,seq,len}"), emit: kma_db  (optional)
        path  "versions.yml"                               , emit: versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTA_CLEAN {

    tag "${meta.id}"

    label 'process_high_memory'

    // conda "bioconda::seqkit=2.8.0 bioconda::kma=1.4.14 conda-forge::pigz=2.8"
    // container "${ workflow.containerEngine == 'singularity' ?
    //     'oras://community.wave.seqera.io/library/seqkit_kma:latest' :
    //     'biocontainers/seqkit:2.8.0' }"

    input:
    tuple val(meta), path(raw_fasta)

    output:
    tuple val(meta), path("*.ccmetagen.fasta.gz")           , emit: cleaned_fasta
    tuple val(meta), path("*_ccmetagen.{name,seq,len,comp}"), emit: kma_db, optional: true
    path  "versions.yml"                                    , emit: versions

    script:
    def args             = task.ext.args             ?: ''
    def min_len          = task.ext.min_len          ?: '300'
    def mask             = task.ext.mask_low_complex ?: '0'
    def dedup            = task.ext.dedup_sequences  ?: '0'
    def keep_intermediate = task.ext.keep_intermediate ?: '0'

    """
    export MIN_LEN="${min_len}"
    export MASK_LOW_COMPLEX="${mask}"
    export DEDUP_SEQUENCES="${dedup}"
    export COMPRESS_FINAL=1
    export KEEP_INTERMEDIATE="${keep_intermediate}"
    export THREADS="${task.cpus}"

    clean_fasta.sh "${raw_fasta}" clean_out/ ${args}

    # ── stage outputs into work dir ────────────────────────────────────────
    mv clean_out/*.ccmetagen*.fasta.gz ./

    # Stage KMA index files if they were produced
    for ext in name seq len comp; do
        find clean_out/ -name "*_ccmetagen.\${ext}" -exec mv {} ./ \\;
    done

    # ── versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/seqkit v//')
        kma:    \$(kma -v 2>&1 | head -1 | sed 's/KMA-//')
        pigz:   \$(pigz --version 2>&1 | head -1 || echo "N/A")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "stub" | gzip > ${prefix}.ccmetagen.fasta.gz
    touch ${prefix}_ccmetagen.name ${prefix}_ccmetagen.seq ${prefix}_ccmetagen.len
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.8.0
        kma: 1.4.14
        pigz: 2.8
    END_VERSIONS
    """
}

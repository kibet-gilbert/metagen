// =============================================================================
// subworkflows/local/validate_pathogen_detection/main.nf
//
// Re-usable subworkflow: KrakenTools extraction + BLASTn validation.
// Exports all processes individually so the standalone validate_pathogen.nf
// can import them directly — mirrors how ccmetagen.nf is structured.
// =============================================================================

include { KRAKENTOOLS_EXTRACTKRAKENREADS } from '../../../modules/local/krakentools/extractkrakenreads/main'
include { BLAST_BLASTN                   } from '../../../modules/local/blast/blastn/main'
include { VALIDATE_HITS                  } from '../../../modules/local/validate_pathogen/validate_hits/main'
include { AGGREGATE_VALIDATION           } from '../../../modules/local/validate_pathogen/aggregate_validation/main'

// =============================================================================
// params defaults set OUTSIDE the workflow body to script level so they actually set params.*
// =============================================================================
params {
    // ── Required ──────────────────────────────────────────────────────────
    samplesheet          = null          // CSV: sample,fastq_1,fastq_2,kraken_output,kraken_report
    outdir               = './results/validation'

    // ── Pathogen target taxids ────────────────────────────────────────────
    // Mononegavirales  : 11157   Flaviviridae        : 11050
    // Vibrio cholerae  : 666     Salmonella Typhi    : 90370
    // Salmonella Typhim: 90371   E. coli             : 562
    // K. pneumoniae    : 573     A. baumannii        : 470
    // P. aeruginosa    : 287     Neisseria gonorrhoeae: 485
    // Haemophilus inf  : 727     S. aureus           : 1280
    // Bunyavirales     : 1980410
    target_taxids        = "11157,11050,1980410,666,90370,90371,562,573,470,287,485,727,1280"

    // ── BLAST database ────────────────────────────────────────────────────
    blast_db_dir         = null          // path to dir containing BLAST DB
    blast_db_name        = 'core_nt'     // DB name (nt | core_nt | custom)
    blast_db_mode        = 'local'       // local | custom | remote

    // ── KrakenTools options ───────────────────────────────────────────────
    include_children     = true
    include_parents      = false
    max_reads_per_taxid  = 50000
    fastq_output         = false

    // ── BLAST options ─────────────────────────────────────────────────────
    blast_task           = 'megablast'
    blast_evalue         = '1e-20'
    blast_perc_id        = 95
    blast_qcov           = 80
    blast_max_seqs       = 5

    // ── Validation thresholds ─────────────────────────────────────────────
    min_blast_hits       = 3
    min_pident           = 95.0
    min_qcovs            = 80.0
    min_bitscore         = 100.0

    // ── Pipeline options ──────────────────────────────────────────────────
    publish_dir_mode     = 'copy'
    skip_blast           = false
    save_extracted_reads = true

    // ── nf-core/HPC options ───────────────────────────────────────────────
    max_cpus             = 48
    max_memory           = '256.GB'
    max_time             = '72.h'
    help                 = false
}

// =============================================================================
// SUBWORKFLOW — imported by the main metagen pipeline via:
//   include { VALIDATE_PATHOGEN_DETECTION } from './subworkflows/local/validate_pathogen_detection/main'
// =============================================================================
workflow VALIDATE_PATHOGEN_DETECTION {

    take:
    ch_input      // tuple val(meta), path(reads), path(kraken_output), path(kraken_report)
    target_taxids // list or comma-string of NCBI taxids
    blast_db_dir  // path (or [] for remote)
    blast_db_name // string

    main:
    ch_versions = Channel.empty()

    // ── STEP 1: ────────────────────────────────────────────────────────
    // STEP 1: Extract reads matching target taxids from Kraken2 classification

    KRAKENTOOLS_EXTRACTKRAKENREADS(ch_input, target_taxids)
    ch_versions = ch_versions.mix(KRAKENTOOLS_EXTRACTKRAKENREADS.out.versions)

    // Combine R1 and R2 into single FASTA per sample for BLAST
    // (BLAST does not handle P.E. natively — we treat each read as a query)
    ch_combined = KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_r1
        .join(KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_r2, remainder: true)
        .map { meta, r1, r2 ->
            tuple(meta, [r1, r2].findAll { it != null && it.name != 'null' })
        }

    // ── STEP 2: ────────────────────────────────────────────────────────
    // Concatenate paired reads into one FASTA for BLAST input
    CONCAT_EXTRACTED_READS(ch_combined)
    ch_versions = ch_versions.mix(CONCAT_EXTRACTED_READS.out.versions)

    ch_nonempty = CONCAT_EXTRACTED_READS.out.fasta
        .filter { meta, fasta -> fasta.countFasta() > 0 }

    // ── STEP 3: ────────────────────────────────────────────────────────
    // STEP 3: BLAST extracted reads against trusted reference DB
    BLAST_BLASTN(ch_nonempty, blast_db_dir, blast_db_name)
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)

    ch_val_input = KRAKENTOOLS_EXTRACTKRAKENREADS.out.summary
        .join(BLAST_BLASTN.out.summary, remainder: true)
        .map { meta, ks, bs ->
            tuple(meta, ks, (bs && bs.name != 'null') ? bs : [])
        }

    // ── STEP 4: ────────────────────────────────────────────────────────
    // STEP 4: Cross-validate Kraken2 vs BLAST evidence
    VALIDATE_HITS(ch_val_input, target_taxids)
    ch_versions = ch_versions.mix(VALIDATE_HITS.out.versions)

    // ── STEP 5: ────────────────────────────────────────────────────────
    // STEP 5: Aggregate Kraken2 and BLAST hits
    AGGREGATE_VALIDATION(VALIDATE_HITS.out.report.collect { m, t -> t })

    emit:
    extracted_reads   = KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_r1
    extracted_r2      = KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_r2
    extract_summary   = KRAKENTOOLS_EXTRACTKRAKENREADS.out.summary
    blast_tsv         = BLAST_BLASTN.out.tsv
    blast_summary     = BLAST_BLASTN.out.summary
    validation_report = VALIDATE_HITS.out.report
    aggregate_tsv     = AGGREGATE_VALIDATION.out.tsv
    versions          = ch_versions
}

// =============================================================================
// PROCESS: Concatenate R1 + R2 extracted FASTA into one query file for BLAST
// =============================================================================
process CONCAT_EXTRACTED_READS {
    tag "${meta.id}"
    label 'process_low'

    conda     "conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/pigz:2.8' :
        'quay.io/biocontainers/pigz:2.8' }"

    publishDir(
        path:    "${params.outdir}/extracted_reads/combined",
        mode:    params.publish_dir_mode ?: 'copy',
        pattern: "*.combined.fasta",
        enabled: params.save_extracted_reads ?: true
    )

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("${meta.id}.extracted.combined.fasta"), emit: fasta
    path  "versions.yml",                                          emit: versions

    script:
    """
    set -euo pipefail
    : > ${meta.id}.extracted.combined.fasta
    i=1
    for fa in ${fastas.join(' ')}; do
        if [[ "\$fa" == *.gz ]]; then
            zcat "\$fa" | awk -v i="\$i" '/^>/{sub(/^>/,">"i"/"); print; next}{print}' \
              >> ${meta.id}.extracted.combined.fasta
        else
            awk -v i="\$i" '/^>/{sub(/^>/,">"i"/"); print; next}{print}' "\$fa" \
              >> ${meta.id}.extracted.combined.fasta
        fi
        i=\$((i+1))
    done
    n_seqs=\$(grep -c '^>' ${meta.id}.extracted.combined.fasta 2>/dev/null || echo 0)
    echo "[INFO] ${meta.id}: combined FASTA has \${n_seqs} sequences"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //')
    END_VERSIONS
    """

    stub:
    """
    printf ">stub\\nACGT\\n" > ${meta.id}.extracted.combined.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: 2.8
    END_VERSIONS
    """
}


// =============================================================================
// Standalone Nextflow pipeline for orthogonal validation of Kraken2-detected
// pathogens using KrakenTools read extraction and BLASTn confirmation.
//
// Designed to run independently from the main metagen pipeline, taking
// pre-existing Kraken2 outputs as inputs via a samplesheet.
//
// Usage:
//   nextflow run subworkflows/local/validate_pathogen_detection/main.nf \
//     -c conf/nextflow.validate_pathogen.config \
//     --samplesheet samplesheet_validation.csv \
//     --outdir results/validation \
//     --target_taxids "562,1639,666,90370,1280,287,573,345073" \
//     --blast_db_dir /export/data/ncbi/blast/db/v5 \
//     --blast_db_name core_nt \
//     -profile ilri_hpc \
//     -resume
// =============================================================================
// ── Parameter defaults ───────────────────────────────────────────────────────
workflow VALIDATE_PATHOGEN_STANDALONE {

    main:
    // ── Print help ───────────────────────────────────────────────────────────
    if (params.help) {
        log.info """
        ╔══════════════════════════════════════════════════════════════════════╗
        ║         Pathogen Detection Validation Pipeline                       ║
        ║         KrakenTools + BLASTn   ·   ILRI Genomics Platform           ║
        ╚══════════════════════════════════════════════════════════════════════╝
    
        USAGE:
          nextflow run validate_pathogen.nf [options]
    
        REQUIRED:
          --samplesheet     CSV file with columns:
                            sample, fastq_1, fastq_2, kraken_output, kraken_report
          --blast_db_dir    Path to directory containing BLAST database
                            (not required if --blast_db_mode remote)
    
        OPTIONAL:
          --outdir          Output directory  [default: ${params.outdir}]
          --target_taxids   Comma-separated NCBI taxids  [default: ${params.target_taxids}]
          --blast_db_name   BLAST DB name  [default: ${params.blast_db_name}]
          --blast_db_mode   local | custom | remote  [default: ${params.blast_db_mode}]
          --include_children  Include descendant taxa  [default: ${params.include_children}]
          --max_reads_per_taxid  Max reads extracted per taxid  [default: ${params.max_reads_per_taxid}]
          --skip_blast      Extract reads only, skip BLAST  [default: ${params.skip_blast}]
          --blast_perc_id   Min BLAST % identity  [default: ${params.blast_perc_id}]
          --blast_evalue    BLAST e-value cutoff  [default: ${params.blast_evalue}]
          --min_blast_hits  Min hits to CONFIRM detection  [default: ${params.min_blast_hits}]
          -resume           Resume a previous run
    
        PROFILES:
          -profile ilri_hpc     SLURM + Apptainer (ILRI HPC)
          -profile singularity  Generic Singularity
          -profile conda        Conda environments
    
        EXAMPLE:
          nextflow run validate_pathogen.nf \\
            -c conf/nextflow.validate_pathogen.config \\
            --samplesheet data/samplesheet_validation.csv \\
            --outdir results/validation/run01 \\
            --target_taxids "666,562,573,470" \\
            --blast_db_dir /export/data/ncbi/blast/db/v5 \\
            --blast_db_name core_nt \\
            -profile ilri_hpc -resume
        """.stripIndent()
        System.exit(0)
    }
    
    // ── Validate required params ─────────────────────────────────────────────
    if (!params.samplesheet) {
        error """
        [VALIDATE_PATHOGEN_STANDALONE] --samplesheet is required.
        Run with --help for usage.
        """
    }
    if (params.blast_db_mode != 'remote' && !params.blast_db_dir) {
        error("ERROR: --blast_db_dir is required unless --blast_db_mode remote.")
    }
    
    // ── Pipeline log header ──────────────────────────────────────────────────
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║         Pathogen Detection Validation Pipeline                       ║
    ║         ILRI Genomics Platform                                       ║
    ╠══════════════════════════════════════════════════════════════════════╣
    ║  Samplesheet    : ${params.samplesheet}
    ║  Outdir         : ${params.outdir}
    ║  Target taxids  : ${params.target_taxids}
    ║  BLAST DB       : ${params.blast_db_dir}/${params.blast_db_name} (${params.blast_db_mode})
    ║  Inc. children  : ${params.include_children}
    ║  Max reads/taxid: ${params.max_reads_per_taxid}
    ║  Skip BLAST     : ${params.skip_blast}
    ╚══════════════════════════════════════════════════════════════════════╝
    """.stripIndent()
    
    // ── Parse target taxids into a list ──────────────────────────────────────
    def taxid_list = params.target_taxids
        .toString()
        .split(',')
        .collect { it.trim() }
        .findAll { it }
    
    // ── Read and validate samplesheet ────────────────────────────────
    ch_input = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            // Validate required columns
            def required_cols = ['sample','fastq_1','fastq_2','kraken_output','kraken_report']
            required_cols.each { col ->
                if (!row.containsKey(col) || !row[col]) {
                    error("Samplesheet row missing required column '${col}': ${row}")
                }
            }

            def meta = [
                id          : row.sample,
                single_end  : false,
                county      : row.county      ?: 'NA',
                site        : row.site        ?: 'NA',
                week        : row.year_week   ?: row.week ?: 'NA'
            ]

            def reads        = [file(row.fastq_1, checkIfExists: true),
                                file(row.fastq_2, checkIfExists: true)]
            def kraken_out   = file(row.kraken_output,  checkIfExists: true)
            def kraken_rep   = file(row.kraken_report,  checkIfExists: true)

            tuple(meta, reads, kraken_out, kraken_rep)
        }

    ch_versions = Channel.empty()

    // ── Run the subworkflow ──────────────────────────────────────────────────
    VALIDATE_PATHOGEN_DETECTION(
        ch_input,
        taxid_list,
        params.blast_db_dir,
        params.blast_db_name
    )
}

// =============================================================================
// Add bare entry workflow required by Nextflow 26.04 strict syntax.
// The -entry flag is no longer supported; a nameless workflow{} block is the
// only supported entry point when running a file directly.
// =============================================================================
def completionSummary() {
    def status = workflow.success ? 'COMPLETED' : 'FAILED'
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║  Pipeline ${status}
    ║  Duration : ${workflow.duration}
    ║  Exit code: ${workflow.exitStatus}
    ║  Outdir   : ${params.outdir}
    ╚══════════════════════════════════════════════════════════╝
    """.stripIndent()
}

workflow {
    VALIDATE_PATHOGEN_STANDALONE()

    // ── On completion ────────────────────────────────────────────────────────
    // Hooks inside workflow{} — one-liners, no nested GStrings
    workflow.onComplete { completionSummary() }
    workflow.onError { err ->
        log.error "Pipeline FAILED: ${err?.message ?: 'see .nextflow.log for details'}"
    }
}



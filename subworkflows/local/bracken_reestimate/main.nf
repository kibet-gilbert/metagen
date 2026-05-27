/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BRACKEN REESTIMATION SUBWORKFLOW
    nf-core/metagenomics pipeline — local subworkflow
    Compatible with Nextflow >= 25.04

    Takes:  Kraken2 kreport files + a Kraken2/Bracken database
    Does:   (1) Optionally builds Bracken DB files if not already present
            (2) Runs Bracken re-estimation on each kreport
            (3) Combines all per-sample Bracken reports into one merged table

    Usage (standalone):
        nextflow run subworkflows/local/bracken_reestimate/main.nf \
            -entry BRACKEN_REESTIMATE_STANDALONE \
            --kreports  "kreports/kraken2/*.kraken2_report.txt" \
            --db        "/path/to/kraken2_bracken_db" \
            --outdir    "results/bracken" \
            --readlen   150 \
            --level     S \
            --threshold 10 \
            -profile singularity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Module imports  (paths relative to the project root)
// ---------------------------------------------------------------------------
include { BRACKEN_BUILD              } from '../../../modules/nf-core/bracken/build/main'
include { BRACKEN_BRACKEN            } from '../../../modules/nf-core/bracken/bracken/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS } from '../../../modules/nf-core/bracken/combinebrackenoutputs/main'

// ---------------------------------------------------------------------------
// SUBWORKFLOW: BRACKEN_REESTIMATE
//
// Params consumed (via ext.args in config, or passed directly):
//   params.bracken_readlen   : read length used in Kraken2 run (default 150)
//   params.bracken_level     : taxonomic level S|G|F|O|C|P (default 'S')
//   params.bracken_threshold : minimum reads required (default 10)
//   params.bracken_build     : whether to run bracken/build step (default false)
//
// Channels in:
//   kreports  : Channel of [ meta, path(kreport) ]
//   db        : path  — single Kraken2/Bracken database directory
// ---------------------------------------------------------------------------
workflow BRACKEN_REESTIMATE {

    take:
    kreports  // channel: [ val(meta), path(kreport) ]
    db        // value:   path  (single db directory, shared across all samples)

    main:

    ch_versions = Channel.empty()

    // ── STEP 1 (optional): build Bracken database files ──────────────────────
    // Only needed when the *.kmer_distrib file is missing from the db directory.
    // Controlled by params.bracken_build (default: false — assume db is ready).

    if ( params.bracken_build ) {
        // bracken/build takes tuple val(meta), path(db)
        ch_db_input = Channel.value( [ [ id: 'bracken_db' ], file(db) ] )

        BRACKEN_BUILD ( ch_db_input )
        ch_versions = ch_versions.mix( BRACKEN_BUILD.out.versions_bracken )

        // After build, the same db directory now contains the kmer_distrib file
        ch_db_ready = BRACKEN_BUILD.out.db.map { _meta, db_path -> db_path }
    } else {
        ch_db_ready = Channel.value( file(db) )
    }

    // ── STEP 2: Bracken re-estimation (per sample) ───────────────────────────
    // bracken/bracken takes: tuple val(meta), path(report)  and  path(db)
    BRACKEN_BRACKEN (
        kreports,
        ch_db_ready
    )
    ch_versions = ch_versions.mix( BRACKEN_BRACKEN.out.versions_bracken )

    // ── STEP 3: Combine all per-sample Bracken reports ────────────────────────
    // bracken/combinebrackenoutputs takes: tuple val(meta), path(reports...)
    // Collect all .txt outputs into a single channel item with a combined meta.

    ch_reports_collected = BRACKEN_BRACKEN.out.reports
        .map { _meta, report -> report }
        .collect()
        .map { reports ->
            [ [ id: 'all_samples' ], reports ]
        }

    BRACKEN_COMBINEBRACKENOUTPUTS ( ch_reports_collected )
    ch_versions = ch_versions.mix( BRACKEN_COMBINEBRACKENOUTPUTS.out.versions_combine_bracken_outputs )

    emit:
    bracken_reports  = BRACKEN_BRACKEN.out.reports         // channel: [ meta, tsv ]  per-sample
    bracken_txt      = BRACKEN_BRACKEN.out.txt             // channel: [ meta, txt ]  per-sample (kreport-style)
    combined         = BRACKEN_COMBINEBRACKENOUTPUTS.out.txt  // channel: [ meta, txt ]  merged table
    versions         = ch_versions
}


// ---------------------------------------------------------------------------
// STANDALONE ENTRY POINT
// Run this subworkflow directly without a parent pipeline:
//
//   nextflow run subworkflows/local/bracken_reestimate/main.nf \
//       -entry BRACKEN_REESTIMATE_STANDALONE  ...
// ---------------------------------------------------------------------------
workflow BRACKEN_REESTIMATE_STANDALONE {

    // ── Parameter defaults (override on the command line or in nextflow.config)
    params.kreports        = null          // REQUIRED: glob to kreport files
    params.db              = null          // REQUIRED: path to kraken2/bracken DB
    params.outdir          = 'results/bracken'
    params.bracken_readlen = 150
    params.bracken_level   = 'S'
    params.bracken_threshold = 10
    params.bracken_build   = false         // set true if kmer_distrib is missing

    // ── Validate required inputs ─────────────────────────────────────────────
    if ( !params.kreports ) {
        error "[BRACKEN_REESTIMATE] --kreports is required. Example: --kreports 'kreports/kraken2/*.kraken2_report.txt'"
    }
    if ( !params.db ) {
        error "[BRACKEN_REESTIMATE] --db is required. Provide path to your Kraken2/Bracken database directory."
    }

    // ── Build kreports channel from glob ─────────────────────────────────────
    // Each report gets a meta map with id derived from the filename.
    // Filename pattern: <SAMPLE_ID>.kraken2_report.txt
    ch_kreports = Channel
        .fromPath( params.kreports, checkIfExists: true )
        .map { report ->
            def sample_id = report.name.replaceAll(/\.kraken2_report\.txt$/, '')
            def meta      = [ id: sample_id ]
            [ meta, report ]
        }

    // ── Run subworkflow ──────────────────────────────────────────────────────
    BRACKEN_REESTIMATE (
        ch_kreports,
        params.db
    )

    // ── Publish outputs ──────────────────────────────────────────────────────
    // BRACKEN_REESTIMATE.out.bracken_reports
    //     .subscribe { meta, tsv ->
    //         log.info "  [bracken] ${meta.id} → ${tsv}"
    //     }

    BRACKEN_REESTIMATE.out.combined
        .subscribe { _meta, combined ->
            log.info "  [combined] → ${combined}"
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subworkflows/local/ccmetagen_taxonomy/main.nf

    CCMetagen taxonomic profiling subworkflow.
    Handles the complete chain from DB preparation through per-sample classification
    and final merged output.

    ┌─────────────────────────────────────────────────────────────────────────────┐
    │  DB preparation (one of four modes — choose exactly one):                   │
    │    MODE 1: blast_export  — export BLAST DB (nt/core_nt) → FASTA_CLEAN      │
    │    MODE 2: build         — KMA-index a pre-existing FASTA via CCMETAGEN_DB  │
    │    MODE 3: download      — wget pre-built archive via CCMETAGEN_DB          │
    │    MODE 4: prebuilt      — use an existing KMA DB directory (skip build)    │
    │                                                                             │
    │  For modes 1–3 set params.ccmetagen_db_mode accordingly.                   │
    │  For mode 4 set  params.ccmetagen_db to the path of the KMA DB directory.  │
    │                                                                             │
    │  Classification chain (all modes):                                          │
    │    KMA_ALIGN  →  CCMETAGEN  →  CCMETAGEN_MERGE                             │
    └─────────────────────────────────────────────────────────────────────────────┘

    Channels in (when used as a subworkflow inside a parent pipeline):
        reads     : Channel of [ val(meta), path(reads) ]
                    meta must contain { id: <sample_id>, single_end: true|false }
        db        : optional Channel value path — pre-built KMA DB dir
                    (ignored when params.ccmetagen_db_mode != 'prebuilt')

    Standalone entry point:
        nextflow run subworkflows/local/ccmetagen_taxonomy/main.nf \\
            -entry CCMETAGEN_TAXONOMY_STANDALONE \\
            -c conf/nextflow.ccmetagen.config \\
            --reads  "fastq/*_{R1,R2}.fastq.gz" \\
            --db     "/data/ccmetagen/RefSeq_kma" \\
            -profile singularity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Module imports  (paths relative to project root)
// ---------------------------------------------------------------------------
include { BLAST_EXPORT  } from '../../../modules/local/blast_export'
include { FASTA_CLEAN   } from '../../../modules/local/fasta_clean'
include { CCMETAGEN_DB  } from '../../../modules/local/ccmetagen/ccmetagen_db/main'
include { KMA_ALIGN     } from '../../../modules/local/kma_align'
include { CCMETAGEN     } from '../../../modules/local/ccmetagen/ccmetagen/main'
include { CCMETAGEN_MERGE   } from '../../../modules/local/ccmetagen/ccmetagen_merge/main'

// ---------------------------------------------------------------------------
// SUBWORKFLOW: CCMETAGEN_TAXONOMY
//
// Parameters consumed (set in config or CLI):
//   params.ccmetagen_db_mode       'prebuilt' | 'blast_export' | 'build' | 'download'
//   params.ccmetagen_db            path to KMA DB dir (prebuilt mode) OR keyword
//                                  (RefSeq | NCBI_nt | NCBI_nt_no_env) for build/download
//   params.ccmetagen_db_action     'buildccmetagen_db' | 'downloadccmetagen_db'
//                                  (only for build/download modes; derived automatically)
//   params.ccmetagen_taxon         taxon for RefSeq build  (default: 'bacteria')
//   params.ccmetagen_keyword       keyword string for download/build modes
//   params.blast_db_prefix         path to BLAST DB prefix  (blast_export mode only)
//   params.blast_db_name           short label for the BLAST DB  (default: 'nt')
//   params.single_end              treat reads as single-end  (default: false)
// ---------------------------------------------------------------------------
workflow CCMETAGEN_TAXONOMY {

    take:
    reads     // channel: [ val(meta), path(reads) ]
    db        // value channel: path(db_dir)  — only used in prebuilt mode

    main:

    ch_versions = Channel.empty()

    // ── STEP 1: DB preparation ───────────────────────────────────────────────

    def db_mode = params.ccmetagen_db_mode ?: 'prebuilt'

    switch ( db_mode ) {

        // ── MODE 1: Export BLAST DB then clean + KMA-index ──────────────────
        case 'blast_export':

            def blast_prefix = params.blast_db_prefix
            if ( !blast_prefix ) {
                error "[CCMETAGEN_TAXONOMY] blast_export mode requires --blast_db_prefix"
            }
            def blast_name = params.blast_db_name ?: 'nt'
            def blast_meta = [ id: blast_name ]

            // Stage the BLAST DB directory as a channel item
            ch_blast_in = Channel.value( [ blast_meta, file(blast_prefix).parent ] )

            BLAST_EXPORT( ch_blast_in )
            ch_versions = ch_versions.mix( BLAST_EXPORT.out.versions )

            FASTA_CLEAN( BLAST_EXPORT.out.raw_fasta )
            ch_versions = ch_versions.mix( FASTA_CLEAN.out.versions )

            // The KMA DB is the directory containing *.name *.seq *.len
            // FASTA_CLEAN emits it via kma_db channel
            ch_db_ready = FASTA_CLEAN.out.kma_db
                .map { meta, db_files ->
                    // db_files may be multiple index files; return their parent dir
                    def db_dir = db_files instanceof List ? db_files[0].parent : db_files.parent
                    [ meta, db_dir ]
                }
            break

        // ── MODE 2: Build KMA DB from keyword (NCBI_nt / RefSeq) ───────────
        case 'build':

            def kw      = params.ccmetagen_keyword ?: params.ccmetagen_db
            def kw_meta = [ id: kw ?: 'ccmetagen_db' ]
            ch_db_build_in = Channel.value( [ kw_meta, file('NO_FILE') ] )

            CCMETAGEN_DB(
                ch_db_build_in,
                'buildccmetagen_db',
                kw
            )
            ch_versions = ch_versions.mix( CCMETAGEN_DB.out.versions )
            ch_db_ready = CCMETAGEN_DB.out.db
            break

        // ── MODE 3: Download pre-built archive ──────────────────────────────
        case 'download':

            def kw      = params.ccmetagen_keyword ?: params.ccmetagen_db
            def kw_meta = [ id: kw ?: 'ccmetagen_db' ]
            ch_db_dl_in = Channel.value( [ kw_meta, file('NO_FILE') ] )

            CCMETAGEN_DB(
                ch_db_dl_in,
                'downloadccmetagen_db',
                kw
            )
            ch_versions = ch_versions.mix( CCMETAGEN_DB.out.versions )
            ch_db_ready = CCMETAGEN_DB.out.db
            break

        // ── MODE 4: Pre-built KMA DB directory (default) ────────────────────
        case 'prebuilt':
        default:

            if ( db ) {
                // passed directly from parent pipeline as a path value channel
                ch_db_ready = db.map { d -> [ [ id: 'prebuilt_db' ], d ] }
            } else if ( params.ccmetagen_db ) {
                ch_db_ready = Channel.value( [ [ id: 'prebuilt_db' ], file(params.ccmetagen_db, checkIfExists: true) ] )
            } else {
                error "[CCMETAGEN_TAXONOMY] prebuilt mode requires --db or --ccmetagen_db pointing to a KMA DB directory"
            }
            break
    }

    // ── STEP 2: KMA alignment (per sample, parallel) ─────────────────────────
    // Combine every sample with the single shared DB
    ch_reads_with_db = reads.combine( ch_db_ready )
        .map { meta, sample_reads, db_meta, db_dir ->
            [ meta, sample_reads, db_meta, db_dir ]
        }

    KMA_ALIGN(
        ch_reads_with_db.map { meta, reads, _db_meta, _db_dir -> [ meta, reads ] },
        ch_reads_with_db.map { _meta, _reads, db_meta, db_dir -> [ db_meta, db_dir ] }
    )
    ch_versions = ch_versions.mix( KMA_ALIGN.out.versions.first() )

    // ── STEP 3: CCMetagen taxonomic classification (per sample) ──────────────
    ch_kma_out = KMA_ALIGN.out.res
        .join( KMA_ALIGN.out.mapstat, by: 0 )

    CCMETAGEN(
        ch_kma_out,
        ch_db_ready
    )
    ch_versions = ch_versions.mix( CCMETAGEN.out.versions.first() )

    // ── STEP 4: Merge all per-sample CSVs into one combined table ────────────
    ch_csvs_collected = CCMETAGEN.out.csv
        .map  { _meta, csv -> csv }
        .collect()
        .map  { csvs -> [ [ id: 'all_samples' ], csvs ] }

    CCMETAGEN_MERGE( ch_csvs_collected )
    ch_versions = ch_versions.mix( CCMETAGEN_MERGE.out.versions )

    emit:
    kma_res          = KMA_ALIGN.out.res              // channel: [ meta, res ]
    kma_mapstat      = KMA_ALIGN.out.mapstat          // channel: [ meta, mapstat ]
    ccmetagen_csv    = CCMETAGEN.out.csv              // channel: [ meta, csv ]
    ccmetagen_dir    = CCMETAGEN.out.results_dir      // channel: [ meta, dir ]
    combined         = CCMETAGEN_MERGE.out.combined   // channel: [ meta, combined_csv ]
    versions         = ch_versions
}

// ---------------------------------------------------------------------------
// STANDALONE ENTRY POINT
// Run directly:
//   nextflow run subworkflows/local/ccmetagen_taxonomy/main.nf \
//       -entry CCMETAGEN_TAXONOMY_STANDALONE \
//       -c conf/nextflow.ccmetagen.config \
//       --reads "fastq/*_{R1,R2}.fastq.gz" \
//       --db    "/data/ccmetagen/RefSeq_kma" \
//       -profile singularity
// ---------------------------------------------------------------------------
workflow CCMETAGEN_TAXONOMY_STANDALONE {

    // ── Defaults ────────────────────────────────────────────────────────────
    params.reads               = null
    params.db                  = null
    params.outdir              = 'results/ccmetagen'
    params.single_end          = false
    params.ccmetagen_db_mode   = 'prebuilt'   // prebuilt | blast_export | build | download
    params.ccmetagen_db        = null
    params.ccmetagen_keyword   = null
    params.ccmetagen_taxon     = 'bacteria'
    params.blast_db_prefix     = null
    params.blast_db_name       = 'nt'

    // ── Validate ────────────────────────────────────────────────────────────
    if ( !params.reads ) {
        error """
        [CCMETAGEN_TAXONOMY_STANDALONE] --reads is required.
        Examples:
          paired-end  : --reads "fastq/*_{R1,R2}.fastq.gz"
          single-end  : --reads "fastq/*.fastq.gz" --single_end
        """
    }

    // ── Build reads channel ─────────────────────────────────────────────────
    if ( params.single_end ) {
        ch_reads = Channel
            .fromPath( params.reads, checkIfExists: true )
            .map { reads ->
                def sample_id = reads.name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                def meta = [ id: sample_id, single_end: true ]
                [ meta, [ reads ] ]
            }
    } else {
        // Paired-end: fromFilePairs groups R1/R2 automatically
        ch_reads = Channel
            .fromFilePairs( params.reads, checkIfExists: true )
            .map { sample_id, reads ->
                def meta = [ id: sample_id, single_end: false ]
                [ meta, reads ]
            }
    }

    // ── DB value channel (prebuilt mode) ────────────────────────────────────
    def db_ch = params.ccmetagen_db
        ? Channel.value( file(params.ccmetagen_db, checkIfExists: true) )
        : Channel.empty()

    // ── Run the subworkflow ──────────────────────────────────────────────────
    CCMETAGEN_TAXONOMY(
        ch_reads,
        db_ch
    )

    // ── Log key outputs ──────────────────────────────────────────────────────
    CCMETAGEN_TAXONOMY.out.ccmetagen_csv
        .subscribe { meta, csv ->
            log.info "  [ccmetagen] ${meta.id} → ${csv}"
        }

    CCMETAGEN_TAXONOMY.out.combined
        .subscribe { _meta, combined ->
            log.info "  [merged]    → ${combined}"
        }
}

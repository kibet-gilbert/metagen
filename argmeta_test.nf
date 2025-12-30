#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Reads:
params.input='samplesheet.csv'
params.single_end=false
// params.reads='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/fastq/COVG0003*con_R{1,2}_001.fastq.gz'
// reads = Channel.fromFilePairs(params.reads, checkIfExists:true)
// println(reads.view())

// hostile indices
params.hostile_db='/export/data/ilri/sarscov2/databases/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401/'
Channel
	.fromPath(params.hostile_db, type: 'dir', checkIfExists:true)
        .map { myHostileDB -> [myHostileDB.getName(), myHostileDB.toString()] }
	.collect()
	.set { hostileDB }
// println(hostileDB.view())

// resfinder database
params.resfinder_db='/export/data/ilri/sarscov2/databases/resfinder_db/'
Channel
	.fromPath(params.resfinder_db, type: 'dir', checkIfExists:true)
	.map { myResFinderDB -> [myResFinderDB.getSimpleName(), myResFinderDB.toString()] }
 	.collect()
	.set { resfinder_db }

// pointfinder database
params.pointfinder_db='/export/data/ilri/sarscov2/databases/pointfinder_db/'
Channel
	.fromPath(params.pointfinder_db, type: 'dir', checkIfExists:true)
	.map { myPointFinderDB -> [myPointFinderDB.getSimpleName(), myPointFinderDB.toString()] }
 	.collect()
	.set { pointfinder_db }

// disinfinder database
params.disinfinder_db='/export/data/ilri/sarscov2/databases/disinfinder_db/'
Channel
	.fromPath(params.disinfinder_db, type: 'dir', checkIfExists:true)
	.map { myDisinFinderDB -> [myDisinFinderDB.getSimpleName(), myDisinFinderDB.toString()] }
 	.collect()
	.set { disinfinder_db }

// RGI
params.carddbcversion = '4.0.1'
params.carddbwversion = '4.0.2'
params.CARD_db='/export/data/ilri/sarscov2/databases/card/localDB'
channel
	.fromPath(params.CARD_db, type: 'dir', checkIfExists:true)
        .map { myCARDDB -> [myCARDDB.getSimpleName(), myCARDDB.toString()] }
	.collect()
	.set { CARD_db }

// MGS2AMR
params.mgs2amr_db='/export/data/ilri/sarscov2/databases/mgs2amr/mgs2amr_assets'
Channel
        .fromPath(params.mgs2amr_db, type: 'dir', checkIfExists: true)
        .map { MGS2AMRDB -> [MGS2AMRDB.getSimpleName(), MGS2AMRDB.toString()] }
        .collect()
        .set { mgs2amr_db }
// MGS2AMR params
params.mgs2amr_step = 4
params.mgs2amr_verbose = 1
params.mgs2amr_compress = true
params.mgs2amr_force = true

// BLASTDB
params.blast_db='/export/data/bio/ncbi/blast/db/v5/'
Channel
        .fromPath(params.blast_db, type: 'dir', checkIfExists: true)
        .map { BLASTDB -> [BLASTDB.getSimpleName(), BLASTDB.toString()] }
        .collect()
        .set { blast_db }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { fromSamplesheet } from 'plugin/nf-validation'
include { FASTQC } from './modules/local/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED} from './modules/local/fastqc.nf'
include { FASTP } from './modules/local/fastp.nf'
include { HOSTILE } from './modules/local/hostile.nf'
include { RESFINDER } from './modules/local/resfinder.nf'
include { CARDRGI } from './modules/local/cardrgi.nf'
include { MGS2AMR } from './modules/local/mgs2amr.nf'

// Validate channels from input samplesheet
def validateInputSamplesheet(sample_id, mode, sr1, sr2, lr) {
    // println "Validating sample_id:${sample_id}, mode:${mode}, SR1:${sr1}, SR2:${sr2}, LR:${lr}"
    if ( !sr2 && !params.single_end ) {
        error("ERROR: Single-end data must be executed with `--single_end`. Note that it is not possible to mix single- and paired-end data in one run! Check input TSV for sample: ${sample_id}")
    }
    if ( sr2 && params.single_end ) {
        error("[nf-core/mag] ERROR: Paired-end data must be executed without `--single_end`. Note that it is not possible to mix single- and paired-end data in one run! Check input TSV for sample: ${sample_id}")
        }
    return [sample_id, mode, sr1, sr2, lr]
}

// Create channels from input file provided through params.input 
// Validate FASTQ input
// println "Input file path: ${params.input}"
ch_samplesheet = Channel
    .fromPath(params.input, type:'file', checkIfExists:true)
    .splitCsv( header: ['sample_id', 'mode', 'sr1', 'sr2', 'lr'], skip: 1 )
    .map { row ->
        // println "Row content: ${row}"
        validateInputSamplesheet(row.sample_id, row.mode, row.sr1, row.sr2, row.lr) 
    }
    .flatMap { sampleEntry ->
        def (sample_id, mode, sr1, sr2, lr) = sampleEntry 
        def files = [sr1, sr2, lr].findAll { it != null && it != '' } // Only keep non-null values
        [[sample_id, mode, sr1, sr2, lr]] 
    }

// Prepare FASTQs channel and separate short and long reads and prepare
ch_raw_short_reads = ch_samplesheet
    .map { sample_id, mode, sr1, sr2, lr ->
                // sample_id.run = sample_id.run == null ? "0" : sample_id.run
                // sample_id.single_end = params.single_end

                if (params.single_end) {
                    return [ sample_id, [ sr1 ] ]
                } else {
                    return [ sample_id, [ sr1, sr2 ] ]
                }
        }
// println(ch_raw_short_reads.view())
ch_raw_long_reads = ch_samplesheet
    .map { sample_id, mode, sr1, sr2, lr ->
                if (lr) {
                    // sample_id.run = sample_id.run == null ? "0" : sample_id.run
                    return [ sample_id, lr ]
                }
         }

// Workflow definition
workflow {
    // ch_samplesheet.view()
    // ch_raw_short_reads.view()
    // FASTQC(ch_raw_short_reads)
    FASTP(ch_raw_short_reads)
    // println(FASTP.out.trimmed_reads.view())
    // FASTQC_TRIMMED(FASTP.out.trimmed_reads)
    HOSTILE(FASTP.out.trimmed_reads,
        hostileDB)
    RESFINDER(HOSTILE.out.reads,
        resfinder_db,
        pointfinder_db,
        disinfinder_db)
    CARDRGI(HOSTILE.out.reads, 
        CARD_db)
    MGS2AMR(HOSTILE.out.reads, 
        params.mgs2amr_step,
        params.mgs2amr_verbose,
        params.mgs2amr_compress,
        params.mgs2amr_force,
        mgs2amr_db)
    // BUILD_KRAKEN2DB(kraken2_human_db)
    // Kraken2Host(FastP.out.trimmed_reads, BUILD_KRAKEN2DB.out.kraken2_db)
    // BUILD_CENTRIFUGEDB(centrifuge_db)
    // CENTRIFUGE(Kraken2Host.out.nohost_reads, BUILD_CENTRIFUGEDB.out.centrifuge_db)
    // BUILD_KronaDB(krona_db)
    // KronaTools(CENTRIFUGE.out.krona_report, BUILD_KronaDB.out.kronaDB)
    // BUILD_KRAKEN2DB(kraken2_db)
    // Kraken2Taxonomy(Kraken2Host.out.nohost_reads, BUILD_KRAKEN2DB.out.kraken2_db)
    // KronaTools(Kraken2Taxonomy.out.kraken_report, BUILD_KronaDB.out.kronaDB)
    // MetaSPAdes(Kraken2Host.out.nohost_reads)
    // println(Kraken2Host.out.nohost_reads.collect().view())
    // Kraken2Contigs(MetaSPAdes.out.contigs,kraken2_db_viral)
    // KronaContigs(Kraken2Contigs.out.krona_report,taxonomy)
    // SAM2FASTQ(arg_bam)
    // BBTOOLSREPAIR(SAM2FASTQ.out.samFile)
    // BUILD_KRAKEN2DB(kraken2_db)
    // Kraken2Taxonomy(SAM2FASTQ.out.reads, BUILD_KRAKEN2DB.out.kraken2_db)
    // BUILD_KronaDB(krona_db)
    // KRONAPLOT(Kraken2Taxonomy.out.krona_report, BUILD_KronaDB.out.kronaDB)
    // BUILD_CENTRIFUGEDB(centrifuge_db)
    // CENTRIFUGE(BBTOOLSREPAIR.out.trimmed_reads, BUILD_CENTRIFUGEDB.out.centrifuge_db)
    // KRONAPLOTB(CENTRIFUGE.out.krona_report, BUILD_KronaDB.out.kronaDB)
}

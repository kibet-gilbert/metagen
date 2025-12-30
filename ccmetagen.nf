#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// params:
// Reads:
// params.reads='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/fastq/*con_R{1,2}_001.fastq.gz'
params.kraken2_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/20231129_k2_nt_20231129/'
Channel
        .fromPath(params.kraken2_db, checkIfExists:true)
        .collect()
        .set { kraken2_db }
// Accepts either a path to a file of a HTTP(S)/FTP link to the file,
// params.ccm_db_input → directory / zip / db-name
params.ccm_db_input='/export/data/ilri/sarscov2/databases/ccmetagen/RefSeq_bf.zip'
Channel
    .value(params.ccm_db_input)
    .map { input ->
        def p = file(input)

        if (p.exists()) {
            if (p.isDirectory()) {
                // PRIORITY 2
                return tuple(p.getSimpleName(), p)
            }
            else if (p.name.endsWith(".zip")) {
                // PRIORITY 1
                def base = p.getSimpleName().replaceAll(/\.zip$/, "")
                return tuple(base, p)
            }
        }
        // PRIORITY 3 – keyword
        return tuple(input.toString(), null)
    }
    .set { ch_ccm_db_in }

// Action channel
// params.db_mode      → "buildccmetagen_db" | "downloadccmetagen_db"
params.db_mode='downloadccmetagen_db'
ch_ccm_db_action = Channel.value(params.db_mode)

// Building CCMetagen RefSeq DB for taxon
params.ccmetagen_taxon='All'

include { FASTP }          from './modules/local/fastp.nf'
include { CCMETAGEN_DB }   from './modules/local/ccmetagen_db.nf'
include { KMA }            from './modules/local/kma.nf'
include { CCMETAGEN }      from './modules/local/ccmetagen.nf'

workflow {

    /*
     * ===============================
     *   Input handling
     * ===============================
     */
    Channel
        .fromFilePairs("${params.reads}", checkIfExists: true)
        .map { sid, reads -> tuple(sid, reads) }
        .set { ch_reads }

    /*
     * ===============================
     *   DB preparation
     * ===============================
     *
     * params.ccm_db_input → directory / zip / db-name
     * params.db_mode      → "buildccmetagen_db" | "downloadccmetagen_db"
     *
     */
    // ch_ccm_db_in     = Channel.value(params.ccm_db_input)
    ch_ccm_db_action = Channel.value(params.db_mode)
    CCMETAGEN_DB(ch_ccm_db_in, ch_ccm_db_action)

    // Output channel
    CCMETAGEN_DB.out.set { ch_ccm_db }


    /*
     * ===============================
     *   FASTP QC
     * ===============================
     */
    FASTP(ch_reads)
    // FASTP.out.set { ch_filtered_reads }


    /*
     * ===============================
     *   KMA alignment
     * ===============================
     */
    KMA( FASTP.out.trimmed_reads, ch_ccm_db )
    KMA.out.set { ch_kma_results }


    /*
     * ===============================
     *   CCMETAGEN profiling
     * ===============================
     */
    CCMETAGEN(ch_kma_results, ch_ccm_db )

}


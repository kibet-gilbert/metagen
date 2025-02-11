#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Reads:
// params.reads='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/fastq/*con_R{1,2}_001.fastq.gz'
params.reads='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/fastq/COVG0003*con_R{1,2}_001.fastq.gz'
reads = Channel.fromFilePairs(params.reads, checkIfExists:true)
// println(reads.view())

// BAM files
params.arg_bam='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/results_rgi/rgi/*.sorted.length_100.bam'
Channel
	.fromPath(params.arg_bam, type:'file', checkIfExists:true)
	.map { myBamFile -> [myBamFile.getSimpleName(), myBamFile] }
	.set { arg_bam }
// arg_bam.view()

// Databases
// RGI
// RGI_db = file('/export/data/ilri/sarscov2/databases/card/localDB', type: 'dir')
params.RGI_db='/export/data/ilri/sarscov2/databases/card/localDB'
channel
	.fromPath(params.RGI_db, type: 'dir', checkIfExists:true)
	.collect()
	.set { RGI_db }
// kraken2
// params.kraken2_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/20231129_k2_nt_20231129.tar.gz'
// params.kraken2_db='https://genome-idx.s3.amazonaws.com/kraken/k2_nt_20231129.tar.gz'
params.kraken2_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/20231129_k2_nt_20231129/'
Channel
	.fromPath(params.kraken2_db, checkIfExists:true)
	.collect()
	.set { kraken2_db }
// Accepts either a path to a file of a HTTP(S)/FTP link to the file
// centrifuge
// params.centrifuge_db='https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz'
params.centrifuge_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/20240222_centrifugeDb_hbpvfa.tar.gz'
Channel
	.fromPath(params.centrifuge_db, checkIfExists:true)
	.map { myCentrifugeDB -> [myCentrifugeDB.getSimpleName(), myCentrifugeDB] }
	.collect()
	.set { centrifuge_db }
// kraken2 human
params.kraken2_human_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/kraken2_human_db.tar.gz'
Channel
	.fromPath(params.kraken2_human_db, checkIfExists:true)
	.collect()
	.set { kraken2_human_db }
// krona Taxonomy
// params.krona_db='https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'
params.krona_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/krona/'
Channel
	.fromPath(params.krona_db, type: 'dir', checkIfExists:true)
	.collect()
	.set { krona_db }

include { FASTQC } from '../modules/local/fastqc.nf'
include { FASTP } from '../modules/local/fastp.nf'
include { KRAKEN2DB as KRAKEN2DBHOST } from '../modules/local/kraken2_db.nf'
include { KRAKEN2FILTER } from '../modules/local/kraken2filter.nf'
include { CARDRGI } from '../modules/local/cardrgi.nf'
include { SAM2FASTQ } from '../modules/local/sam2fastq.nf'
include { BBTOOLSREPAIR } from '../modules/local/bbtools_repair.nf'
include { CENTRIFUGEDB } from '../modules/local/centrifuge_db.nf'
include { CENTRIFUGE } from '../modules/local/centrifuge.nf'
include { KRAKEN2DB } from '../modules/local/kraken2_db.nf'
include { KRAKEN2} from '../modules/local/kraken2.nf'
include { KRONADB } from '../modules/local/krona_db.nf'
include { KRONA as KRONA_KRAKEN2 } from '../modules/local/krona.nf'
include { KRONA as KRONA_CENTRIFUGE } from '../modules/local/krona.nf'
include { KRAKEN2FASTA } from '../modules/local/kraken2fasta.nf'
include { KRONA as KRONA_CONTIGS } from '../modules/local/krona.nf'
include { METASPADES } from '../modules/local/metaspades.nf'

// '../modules/local/bowtie2_db.nf'
// '../modules/local/bowtie2.nf'
// '../modules/local/megahit.nf'
// '../modules/local/metaspades.nf'


// Workflow definition
workflow {
    FASTQC(reads)
    FASTP(reads)
    // KRAKEN2DBHOST(kraken2_human_db)
    // KRAKEN2FILTER(FASTP.out.trimmed_reads, KRAKEN2DBHOST.out.kraken2_db)
    CARDRGI(FASTP.out.trimmed_reads, RGI_db)
    // SAM2FASTQ(arg_bam)
    // BBTOOLSREPAIR(SAM2FASTQ.out.samFile)
    // CENTRIFUGEDB(centrifuge_db)
    // CENTRIFUGE(BBTOOLSREPAIR.out.trimmed_reads, CENTRIFUGEDB.out.centrifuge_db)
    // KRONADB(krona_db)
    // KRONA_CENTRIFUGE(CENTRIFUGE.out.krona_report, KRONADB.out.kronaDB)
    // KRAKEN2DB(kraken2_db)
    // KRAKEN2(BBTOOLSREPAIR.out.trimmed_reads, KRAKEN2DB.out.kraken2_db)
    // KRONA_KRAKEN2(KRAKEN2.out.krona_report, KRONADB.out.kronaDB)
    // METASPADES(KRAKEN2FILTER.out.nohost_reads)
    // println(KRAKEN2FILTER.out.nohost_reads.view())
    // println(KRAKEN2FILTER.out.nohost_reads.collect().view())
    // KRAKEN2FASTA(METASPADES.out.contigs, KRAKEN2DB.out.kraken2_db)
    // KRONA_CONTIGS(KRAKEN2FASTA.out.krona_report, KRONADB.out.kronaDB)
}


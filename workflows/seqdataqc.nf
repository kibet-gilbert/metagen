#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// params:
// PhiX
// params.phix='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz'
params.phix='/export/data/ilri/sarscov2/databases/metagenomicsDBs/GCA_002596845.1_ASM259684v1_genomic.fna'
Channel
	.fromPath(params.phix, checkIfExists:true)
	.map { PhiX -> [PhiX.getSimpleName(), PhiX] }
	.collect()
	.set { phix }
// println(phix.view())
params.host_removal_save_ids = true
// params.megahit_fix_cpu_1 = false
// Reads:
// params.reads='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/fastq/*con_R{1,2}_001.fastq.gz'
reads = Channel.fromFilePairs(params.reads, checkIfExists:true)
// println(reads.view())

// BAM files
// params.arg_bam='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/results_rgi/rgi/*.sorted.length_100.bam'
// Channel
// 	.fromPath(params.arg_bam, type:'file', checkIfExists:true)
// 	.map { myBamFile -> [myBamFile.getSimpleName(), myBamFile] }
// 	.set { arg_bam }
// arg_bam.view()

// Databases
// RGI
// RGI_db = file('/export/data/ilri/sarscov2/databases/card/localDB', type: 'dir')
// params.RGI_db='/export/data/ilri/sarscov2/databases/card/localDB'
// channel
// 	.fromPath(params.RGI_db, type: 'dir', checkIfExists:true)
// 	.collect()
// 	.set { RGI_db }
// kraken2
// params.kraken2_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/20231129_k2_nt_20231129.tar.gz'
// params.kraken2_db='https://genome-idx.s3.amazonaws.com/kraken/k2_nt_20231129.tar.gz'
params.kraken2_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/k2_core_nt_20241228/'
Channel
	.fromPath(params.kraken2_db, checkIfExists:true)
	.collect()
	.set { kraken2_db }
// Accepts either a path to a file of a HTTP(S)/FTP link to the file
// centrifuge
// params.centrifuge_db='https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz'
// params.centrifuge_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/20240222_centrifugeDb_hbpvfa.tar.gz'
// Channel
// 	.fromPath(params.centrifuge_db, checkIfExists:true)
// 	.map { myCentrifugeDB -> [myCentrifugeDB.getSimpleName(), myCentrifugeDB] }
// 	.collect()
// 	.set { centrifuge_db }
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	IMPORT MODULES / LOCAL / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BOWTIE2_DB } from '../modules/local/bowtie2_db.nf'
include { BOWTIE2 as BOWTIE2 } from '../modules/local/bowtie2.nf'
include { FASTQC as FASTQC } from '../modules/local/fastqc.nf'
include { FASTP as FASTP } from '../modules/local/fastp.nf'
include { KRAKEN2DB } from '../modules/local/kraken2_db.nf'
include { KRAKEN2 } from '../modules/local/kraken2.nf'
include { KRONADB } from '../modules/local/krona_db.nf'
include { KRONA } from '../modules/local/krona.nf'

// Workflow definition
workflow {
    //BOWTIE2_DB(phix)
    //BOWTIE2(reads, BOWTIE2_DB.out.index)
    FASTQC(reads)
    FASTP(reads)
    KRAKEN2DB(kraken2_db)
    KRAKEN2(FASTP.out.trimmed_reads, KRAKEN2DB.out.kraken2_db)
    KRONADB(krona_db)
    KRONA(KRAKEN2.out.krona_report, KRONADB.out.kronaDB)
    //println(FASTP.out.reads.collect().view())
}


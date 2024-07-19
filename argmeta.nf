#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Reads:
// params.reads='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/results_rgi/rgi/fastq/*_R{1,2}.fastq.gz'
// params.reads='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/fastq/*con_R{1,2}_001.fastq.gz'
// params.reads='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/fastq/*con_R{1,2}_001.fastq.gz'
params.reads='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/fastq/COVG0003*con_R{1,2}_001.fastq.gz'
reads = Channel.fromFilePairs(params.reads, checkIfExists:true)
// println(reads.view())

// BAM files
//params.arg_bam='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/results_rgi/rgi/*_S0.sorted.length_100.bam'
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

process FastQC {
    // publishDir "${params.outputDir}/fastqc", mode: 'copy'
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip

    script:
    """
    fastqc -t $task.cpus -o . ${reads}
    """
}

process FastP {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*-trim*"), emit: trimmed_reads
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.json"), emit: json
    tuple val(sample_id), path("*.log"), emit: log


    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${sample_id}-trim_R1.fastq.gz \\
        --out2 ${sample_id}-trim_R2.fastq.gz \\
        --html ${sample_id}.fastp.html \\
	--json ${sample_id}.fastp.json \\
        --failed_out ${sample_id}_failed.fastq.gz \\
        --thread $task.cpus \\
        --detect_adapter_for_pe \\
        --dedup \\
        |& tee -a -i ${sample_id}.fastp.log
    """
}

process Kraken2Host {
    // publishDir "${params.outputDir}/kraken2", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(trimmed_reads)
    path("database/*")

    output:
    tuple val(sample_id), path("*.nohost{.,_}*"), emit: nohost_reads
    tuple val(sample_id), path("*.report.txt"), emit: kraken_hostmap_report1
    tuple val(sample_id), path("*.kraken2.out"), emit: kraken_hostmap_report2

    script:
    """
    kraken2 --db database/ --threads $task.cpus \\
            --unclassified-out ${sample_id}.nohost#.fastq \\
            --classified-out ${sample_id}.host#.fastq \\
            --report ${sample_id}.kraken2.report.txt \\
            --output ${sample_id}.kraken2.out \\
            --gzip-compressed --report-zero-counts \\
            --paired $trimmed_reads
    """
}

process CARDRGI {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    path(RGI_db)

    output:
    path("*allele_mapping_data*")
    path("*gene_mapping_data*")
    path("*.temp.sam.temp.*")
    path("*.sorted.length_100.*")
    path("*overall_mapping_stats*")
    path("*.sorted.length_100*")

    script:
    """
    rgi bwt \\
	-1 ${trimmed_reads[0]} \\
	-2 ${trimmed_reads[1]} \\
	-a kma \\
	-n 20 \\
	-o ${sample_id} \\
	--local \\
	--include_other_models \\
	--include_wildcard
    """
}

process SAM2FASTQ {
    // debug true
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(argbamfile)

    output:
    path("*.flagstat.txt"), emit: flagstat_report
    tuple val(sample_id), path("*.fastq.gz"), emit: reads
    tuple val(sample_id), path("*.sam"), emit: samFile
    tuple val(sample_id), path("*.log"), emit: log

    script:
    //myBamFile = file()
    samtoolsSIFRun='apptainer run /home/gkibet/bioinformatics/github/metagenomics/scripts/samtools/samtools_1.20--ad906e74fde1812b.sif samtools'
    """
    ${samtoolsSIFRun} flagstat ${argbamfile} \\
	-@ $task.cpus \\
	|& tee ${sample_id}.flagstat.txt 
    ${samtoolsSIFRun} view ${argbamfile} \\
	-hu \\
	-G 0xd \\
	-@ $task.cpus |
    ${samtoolsSIFRun} sort -n - \\
	-@ $task.cpus \\
	-o ${sample_id}.sortednFiltered.argSeqs.sam \\
	|& tee -a ${sample_id}.samtoolsfastq.log
    ${samtoolsSIFRun} fastq ${sample_id}.sortednFiltered.argSeqs.sam \\
	-@ $task.cpus \\
	-1 ${sample_id}.argSeqs_R1.fastq.gz \\
	-2 ${sample_id}.argSeqs_R2.fastq.gz \\
	-n \\
	|& tee -a ${sample_id}.samtoolsfastq.log
	##-0 ${sample_id}.argSeqs_RO.fastq.gz \\
	##-s ${sample_id}.argSeqs_SR.fastq.gz \\
    """
}

process BBTOOLSREPAIR {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(samFile)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: trimmed_reads

    script:
    """
    repair.sh in=${samFile} \\
	out1=${sample_id}.argSeqs.R1.fastq \\
	out2=${sample_id}.argSeqs.R2.fastq \\
	outs=${sample_id}.argSeqs.SE.fastq \\
	repair=t
    pigz ${sample_id}.argSeqs.R1.fastq
    pigz ${sample_id}.argSeqs.R2.fastq
    pigz ${sample_id}.argSeqs.SE.fastq
    """
}

process BUILD_KRAKEN2DB {
    // debug true
    
    input:
    path(kraken2_db)
    
    output:
    path("database/*.k2d"), emit: kraken2_db

    script:
    if ( kraken2_db.extension in ['gz', 'tgz'] )
	"""
	mkdir -p ./{database,db_tmp}
	tar -xf ${kraken2_db} -C db_tmp
	mv `find db_tmp/ -name "*.k2d"` database/
	"""
    else if ( kraken2_db.isDirectory() ) 
	"""
	mkdir -p database
	cp `find ${kraken2_db}/ -name "*.k2d"` database/
	"""
    else
	error "Path or Link to a kraken2 database not provided in ${params.kraken2_db}"
}

process Kraken2Taxonomy {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    path("database/*")

    output:
    //tuple val(sample_id), path("*.classified{.,_}*"), emit: classified_reads
    tuple val(sample_id), path("*_kreport.txt"), emit: kraken_kreport
    tuple val(sample_id), path("*_kraken2.out"), emit: kraken_report
    tuple val(sample_id), path("*_kraken2.krona"), emit: krona_report

    script:
    """
     kraken2 --db database/ \\
            --threads $task.cpus \\
            --unclassified-out ${sample_id}.unclassified#.fastq \\
            --classified-out ${sample_id}.classified#.fastq \\
            --report ${sample_id}_tax_kreport.txt \\
            --output ${sample_id}_tax_kraken2.out \\
            --report-zero-counts \\
            --gzip-compressed \\
            --paired ${trimmed_reads[0]} ${trimmed_reads[1]}
    cat ${sample_id}_tax_kraken2.out | cut -f 2,3 > ${sample_id}_kraken2.krona
    """
}

process BUILD_CENTRIFUGEDB {
    // debug true
    
    input:
    tuple val(centrifugeDBName), path(centrifuge_db)
    
    output:
    tuple val(centrifugeDBName), path("database/*.cf"), emit: centrifuge_db

    script:
    if ( centrifuge_db.extension in ['gz', 'tgz'] )
	"""
	mkdir -p ./{database,db_tmp}
	tar -xf ${centrifuge_db} -C db_tmp
	mv `find db_tmp/ -name "*.cf"` database/
	"""
    else if ( centrifuge_db.isDirectory() ) 
	"""
	mkdir -p database
	cp `find ${centrifuge_db}/ -name "*.cf"` database/
	"""
    else
	error "Path or Link to a kraken2 database not provided in ${params.centrifuge_db}"
}

process CENTRIFUGE {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    tuple val(centrifugeDBName), path("database/*")

    output:
    tuple val(sample_id), path("*_report.txt"), emit: centrifuge_report
    tuple val(sample_id), path("*_results.txt"), emit: centrifuge_results
    tuple val(sample_id), path("*_kreport2.txt"), emit: kraken_report
    tuple val(sample_id), path("*_centrifuge.krona"), emit: krona_report

    script:
    centrifugeSIFRun='apptainer run /home/gkibet/bioinformatics/github/metagenomics/scripts/apptainer/centrifuge_imp_python.sif'
    """
    #${centrifugeSIFRun} \\
    centrifuge -x database/"${centrifugeDBName}" \\
           -p ${task.cpus} \\
           --report-file ${sample_id}_report.txt \\
           -S ${sample_id}_results.txt \\
           -1 ${trimmed_reads[0]} \\
           -2 ${trimmed_reads[1]}
           #-r ${trimmed_reads[2]}
    #${centrifugeSIFRun} \\
    centrifuge-kreport -x database/"${centrifugeDBName}" \\
           ${sample_id}_results.txt > ${sample_id}.centrifuge_kreport2.txt \\
           2> >(tee ${sample_id}.centrifuge_kreport2_skipped.txt >&2)
    cat ${sample_id}_results.txt | cut -f 1,3 > ${sample_id}_centrifuge.krona
    """
}

process BUILD_KronaDB {
    // debug true
    
    input:
    path(taxonomyDB)
    
    output:
    path("./taxonomy/"), emit: kronaDB

    script:
    if ( taxonomyDB.extension in ['gz', 'tgz'] )
	"""
	mkdir -p ./taxonomy
	tar -xf ${taxonomyDB} -C ./taxonomy/
	ktUpdateTaxonomy.sh --only-build ./taxonomy/
	"""
    else if ( taxonomyDB.isDirectory() ) 
	"""
	mkdir -p ./taxonomy
	cp -rf ${taxonomyDB}/* ./taxonomy/
	"""
    else
	error "Path or Link to a krona database not provided in ${params.krona_db}"
}

process KRONAPLOT {
    tag "$sample_id"  
    
    input:
    tuple val(sample_id), path(krona_report)
    path(taxonomy)
    
    output:
    tuple val(sample_id), path("*.html"), emit: krona_html

    script:
    """    
    ktImportTaxonomy -tax ${taxonomy} \\
                     -o ${sample_id}_Kraken2.krona.html \\
                     ${krona_report}
    """
}

process KRONAPLOTB {
    tag "$sample_id"  
    
    input:
    tuple val(sample_id), path(krona_report)
    path(taxonomy)
    
    output:
    tuple val(sample_id), path("*.html"), emit: krona_html

    script:
    """    
    ktImportTaxonomy -tax ${taxonomy} \\
                     -o ${sample_id}_Centrifuge.krona.html \\
                     ${krona_report}
    """
}

process MetaSPAdes {
    tag "$sample_id"     
    
    input:
    tuple val(sample_id), path(nohost_reads)

    output:
    tuple val(sample_id), path("*_contigs.fasta"), emit: contigs
    tuple val(sample_id), path("*_scaffolds.fasta"), emit: scaffolds
    tuple val(sample_id), path("*_graph.gfa"), emit: graphgfa
    tuple val(sample_id), path("*.log"), emit: spadeslog

    script:
    """
    spades.py \\
        --meta \\
        --threads $task.cpus \\
        --memory 8 \\
        -1 ${nohost_reads[0]} \\
        -2 ${nohost_reads[1]} \\
        -o .
    mv assembly_graph_with_scaffolds.gfa SPAdes-${sample_id}_graph.gfa
    mv scaffolds.fasta SPAdes-${sample_id}_scaffolds.fasta
    mv contigs.fasta SPAdes-${sample_id}_contigs.fasta
    mv spades.log SPAdes-${sample_id}.log
    """
}

process Kraken2Contigs {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(contigs)
    path(kraken2_db_viral)

    output:
    tuple val(sample_id), path("*_kreport.txt"), emit: kraken_viral_report1
    tuple val(sample_id), path("*_kraken2.out"), emit: kraken_viral_report2
    tuple val(sample_id), path("*_kraken2.krona"), emit: krona_report

    script:
    """
    kraken2 --db ${kraken2_db_viral} --threads $task.cpus \\
            --report ${sample_id}_taxContigs_kreport.txt \\
            --output ${sample_id}_taxContigs_kraken2.out \\
            --report-zero-counts $contigs
    cat ${sample_id}_taxContigs_kraken2.out | cut -f 2,3 > ${sample_id}_taxContigs_kraken2.krona
    """
}

process KronaContigs {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(kraken_viral_report2)
    path(taxonomy)

    output:
    tuple val(sample_id), path("*.html"), emit: kronacontigs_html

    script:
    """
    ktImportTaxonomy -tax ${taxonomy} \
                     -o ${sample_id}_taxContigs_taxonomy.krona.html \
                     ${sample_id}_taxContigs_kraken2.krona
    """
}

// include { KRONAPLOT as KRONAPLOT_kraken } from './argmeta.nf'
// include { KRONAPLOT as KRONAPLOT_centrifuge } from './argmeta.nf'

// Workflow definition
workflow {
    // FastP(reads)
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
    // println(Kraken2Host.out.nohost_reads.view())
    // println(Kraken2Host.out.nohost_reads.collect().view())
    // Kraken2Contigs(MetaSPAdes.out.contigs,kraken2_db_viral)
    // KronaContigs(Kraken2Contigs.out.krona_report,taxonomy)
    SAM2FASTQ(arg_bam)
    BBTOOLSREPAIR(SAM2FASTQ.out.samFile)
    BUILD_KRAKEN2DB(kraken2_db)
    Kraken2Taxonomy(SAM2FASTQ.out.reads, BUILD_KRAKEN2DB.out.kraken2_db)
    BUILD_KronaDB(krona_db)
    KRONAPLOT(Kraken2Taxonomy.out.krona_report, BUILD_KronaDB.out.kronaDB)
    BUILD_CENTRIFUGEDB(centrifuge_db)
    CENTRIFUGE(BBTOOLSREPAIR.out.trimmed_reads, BUILD_CENTRIFUGEDB.out.centrifuge_db)
    KRONAPLOTB(CENTRIFUGE.out.krona_report, BUILD_KronaDB.out.kronaDB)
}


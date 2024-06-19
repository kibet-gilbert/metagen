#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// params.reads='/home/gkibet/bioinformatics/metagenomics/data/2022-06-09_run01_nextseq_metagen/results_rgi/rgi/fastq/*_R{1,2}.fastq.gz'
// reads = Channel.fromFilePairs(params.reads, checkIfExists:true)
// println(reads.view())
//params.arg_bam='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/results_rgi/rgi/*_S0.sorted.length_100.bam'
params.arg_bam='/home/gkibet/bioinformatics/github/metagenomics/data/2022-06-09_run01_nextseq_metagen/results_rgi/rgi/*.sorted.length_100.bam'
Channel
	.fromPath(params.arg_bam, type:'file', checkIfExists:true)
	.map { myBamFile -> [myBamFile.getSimpleName(), myBamFile] }
	.set { arg_bam }
// arg_bam.view()
RGI_db = file('/export/data/ilri/sarscov2/databases/card/localDB', type: 'dir')
params.RGI_db='/export/data/ilri/sarscov2/databases/card/localDB'
channel
	.fromPath(params.RGI_db, type: 'dir', checkIfExists:true)
	.set { RGI_db }
// params.kraken2_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/20231129_k2_nt_20231129.tar.gz'
// params.kraken2_db='https://genome-idx.s3.amazonaws.com/kraken/k2_nt_20231129.tar.gz'
params.kraken2_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/20231129_k2_nt_20231129/'
kraken2_db = Channel.fromPath(params.kraken2_db, checkIfExists:true)
// Accepts either a path to a file of a HTTP(S)/FTP link to the file

params.kraken2_human_db='/export/data/ilri/sarscov2/databases/metagenomicsDBs/kraken2_human_db.tar.gz'
kraken2_human_db = Channel.fromPath(params.kraken2_human_db, checkIfExists:true)
params.downloadKronaDB = 'false'
params.krona_db_taxonomy = '/export/data/ilri/sarscov2/databases/metagenomicsDBs/krona/'
krona_db_taxonomy = Channel.fromPath(params.krona_db_taxonomy, type: 'dir', checkIfExists:true)
// params.outputDir = './results'

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
    path(kraken2_db_human)

    output:
    tuple val(sample_id), path("*.nohost{.,_}*"), emit: nohost_reads
    tuple val(sample_id), path("*.report.txt"), emit: kraken_hostmap_report1
    tuple val(sample_id), path("*.kraken2.out"), emit: kraken_hostmap_report2

    script:
    """
    kraken2 --db ${kraken2_db_human} --threads $task.cpus \\
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

    script:
    """
    #KRAKEN2_DB_PATH=./database/
    #apptainer run ~/bioinformatics/github/metagenomics/scripts/kraken2_sing/kraken2_latest.sif 
     kraken2 --db database/ \\
            --threads $task.cpus \\
            --unclassified-out ${sample_id}.unclassified#.fastq \\
            --classified-out ${sample_id}.classified#.fastq \\
            --report ${sample_id}_tax_kreport.txt \\
            --output ${sample_id}_tax_kraken2.out \\
            --report-zero-counts \\
            --gzip-compressed \\
            --paired ${trimmed_reads[0]} ${trimmed_reads[1]}
    """
}

process BUILD_KronaDB {
    // debug true
    
    input:
    path(taxonomyDB)
    
    output:
    path("./taxonomy/"), emit: kronaDB

    script:
    if (params.downloadKronaDB == 'true' )
	"""
	mkdir ./taxonomy/
	wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O ./taxonomy/taxdump.tar.gz
	ktUpdateTaxonomy.sh --only-build ./taxonomy/
	"""
    if (params.downloadKronaDB == 'false' )
	"""
	mkdir ./taxonomy/
	cp -rf ${taxonomyDB}/* ./taxonomy/
	"""
}

process KronaTools {
    tag "$sample_id"  
    
    input:
    tuple val(sample_id), path(kraken_report)
    path(taxonomy)
    
    output:
    tuple val(sample_id), path("*.html"), emit: krona_html

    script:
    """    
    cat $kraken_report | cut -f 2,3 > ${sample_id}_kraken2.krona
    ktImportTaxonomy -tax ${taxonomy} \\
                     -o ${sample_id}_taxonomy.krona.html \\
                     ${sample_id}_kraken2.krona
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

    script:
    """
    kraken2 --db ${kraken2_db_viral} --threads $task.cpus \\
            --report ${sample_id}_taxContigs_kreport.txt \\
            --output ${sample_id}_taxContigs_kraken2.out \\
            --report-zero-counts $contigs
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
    cat $kraken_viral_report2 | cut -f 2,3 > ${sample_id}_taxContigs_kraken2.krona
    ktImportTaxonomy -tax ${taxonomy} \
                     -o ${sample_id}_taxContigs_taxonomy.krona.html \
                     ${sample_id}_taxContigs_kraken2.krona
    """
}

process SAM2FASTQ {
    // debug true
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(argbamfile)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: reads
    tuple val(sample_id), path("*.log"), emit: log

    script:
    //myBamFile = file()
    """
    samtools fastq ${argbamfile} \\
	-@ $task.cpus \\
	-1 ${sample_id}.argSeqs_R1.fastq.gz \\
	-2 ${sample_id}.argSeqs_R2.fastq.gz \\
	-0 ${sample_id}.argSeqs_SE.fastq.gz \\
	-n \\
	|& tee -a -i ${sample_id}.samtoolsfastq.log
    """
}

// Workflow definition
workflow {
    // FastP(reads)
    // Kraken2Host(FastP.out.trimmed_reads, kraken2_db_human)
    // Kraken2Taxonomy(Kraken2Host.out.nohost_reads, kraken2_db_viral)
    // KronaTools(Kraken2Taxonomy.out.kraken_viral_report2, taxonomy)
    // MetaSPAdes(Kraken2Host.out.nohost_reads)
    // println(Kraken2Host.out.nohost_reads.view())
    // println(Kraken2Host.out.nohost_reads.collect().view())
    // Kraken2Contigs(MetaSPAdes.out.contigs,kraken2_db_viral)
    // KronaContigs(Kraken2Contigs.out.kraken_viral_report2,taxonomy)
    SAM2FASTQ(arg_bam)
    BUILD_KRAKEN2DB(params.kraken2_db)
    Kraken2Taxonomy(SAM2FASTQ.out.reads, BUILD_KRAKEN2DB.out.kraken2_db)
    BUILD_KronaDB(params.krona_db_taxonomy)
    KronaTools(Kraken2Taxonomy.out.kraken_report, BUILD_KronaDB.out.kronaDB)
}


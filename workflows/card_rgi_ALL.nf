#!/usr/bin/env nextflow
#SBATCH -p batch
#SBATCH -J cardrgi
#SBATCH --qos normal
#SBATCH -n 16

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

module load fastp/0.22.0
module load fastqc/0.11.9

DBDIR=/export/data/ilri/sarscov2/databases/card/
INDIR=$PWD
src_dir=${INDIR##*/}
WORKDIR=/var/scratch/$USER/${src_dir}
OUTDIR=${WORKDIR}/results_rgi/

mkdir -p ${OUTDIR}/{fastqc,fastp,rgi}
cd ${WORKDIR}
rm -rf ${WORKDIR}/localDB
cp -rf ${DBDIR}localDB ${WORKDIR}/

#Create or setup a samplesheet
source /home/gkibet/bioinformatics/github/metagenomics/scripts/process_files.sh
#FRead_suffix="con_R1_001.fastq.gz"
#RRead_suffix="con_R2_001.fastq.gz"

#samplesheet_gen -f ${FRead_suffix} -r ${RRead_suffix}   


{
	read
	while IFS="" read -r p || [ -n "$p" ];
	do 
		prefix=$(echo $p | cut -f1 -d,);
		fread=$(echo $p | cut -f3 -d,);
		rread=$(echo $p | cut -f4 -d,); 
		echo -e "Processing $prefix\n $fread\n $rread"; #ls $prefix; ls $fread; ls $rread; #done <${INDIR}/fastq/samplesheet.csv
		## Assessing Read Quality using fastqc before quality trimming
		#fastqc -t 4 \
		#-o $OUTDIR/fastqc/ \
		#${fread} \
		#${rread} \
	
		## Quality Trimming fastq files with fastp	
		fastp \
		--in1 ${fread} \
		--in2 ${rread} \
		--out1 $OUTDIR/fastp/${prefix}_1.trim.fastq.gz \
		--out2 $OUTDIR/fastp/${prefix}_2.trim.fastq.gz \
		--json $OUTDIR/fastp/${prefix}.fastp.json \
		--html $OUTDIR/fastp/${prefix}.fastp.html \
		--failed_out $OUTDIR/fastp/${prefix}_fail.fastq.gz \
		--thread 10 \
		--detect_adapter_for_pe \
		--qualified_quality_phred 20 \
		--cut_mean_quality 20 \
		--length_required 15 \
		2> $OUTDIR/fastp/${prefix}.fastp.log
	
		## Assessing Read Quality after quality trimming 	
		#fastqc -t 4 \
		#	-o $OUTDIR/fastqc/ \
		#	$OUTDIR/fastp/${prefix}_1.trim.fastq.gz \
		#	$OUTDIR/fastp/${prefix}_2.trim.fastq.gz \
	
		## RGI ARG gene identification
		apptainer run ~/bioinformatics/github/metagenomics/scripts/card/rgi_6.0.3--pyha8f3691_0.sif rgi bwt \
		-1 $OUTDIR/fastp/${prefix}_1.trim.fastq.gz \
		-2 $OUTDIR/fastp/${prefix}_2.trim.fastq.gz \
		-a kma \
		-n 20 \
		-o $OUTDIR/rgi/${prefix} \
		--local \
		--include_other_models \
		--include_wildcard
	done
}<${INDIR}/samplesheet.csv

## Copy Analysis results to input directory
cp -rf ${OUTDIR} ${INDIR}
if [ $? -eq 0 ]
then
	echo -e "\tCopying successful..."
elif [ $? -ne 0 ]
then
	echo -e "\n\tCopying was NOT successful..."
fi


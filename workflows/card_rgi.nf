process CARD_RGI {
	input:
	output:
	script:
	"""
	## RGI ARG gene identification
	rgi bwt \
	-1  \
	-2 $OUTDIR/fastp/${prefix}_2.trim.fastq.gz \
	-a kma \
	-n 20 \
	-o $OUTDIR/rgi/${prefix} \
	--local \
	--include_other_models \
	--include_wildcard
	"""
}
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


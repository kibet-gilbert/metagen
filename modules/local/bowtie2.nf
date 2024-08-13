/*
 * Bowtie2 for read removal
 */
process BOWTIE2 {
    tag "$sample_id"

    // conda "bioconda::bowtie2=2.4.2"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/bowtie2:2.4.2--py38h1c8e9b9_1' :
    //     'biocontainers/bowtie2:2.4.2--py38h1c8e9b9_1' }"

    input:
    tuple val(sample_id), path(reads)
    tuple val(PhiXAccession), path(index)

    output:
    tuple val(sample_id), path("*.unmapped*.fastq.gz"), emit: reads
    path "*.mapped*.read_ids.txt", optional:true, emit: read_ids
    tuple val(sample_id), path("*.bowtie2.log"), emit: log

    script:
    // def args = task.ext.args ?: ''
    // def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    // def save_ids = (args2.contains('--host_removal_save_ids')) ? "Y" : "N"
    """
    bowtie2 -p ${task.cpus} \
            -x ${index[0].getSimpleName()} \
            -1 "${reads[0]}" -2 "${reads[1]}" \
            $args \
            --un-gz ${prefix}.unmapped.fastq.gz \
            --al-gz ${prefix}.mapped.fastq.gz \
            1> /dev/null \
            2> ${prefix}.bowtie2.log
    gunzip -c ${prefix}.mapped_1.fastq.gz | awk '{if(NR%4==1) print substr(\$0, 2)}' | LC_ALL=C sort > ${prefix}.mapped_1.read_ids.txt
    gunzip -c ${prefix}.mapped_2.fastq.gz | awk '{if(NR%4==1) print substr(\$0, 2)}' | LC_ALL=C sort > ${prefix}.mapped_2.read_ids.txt
    #rm -f ${prefix}.mapped.fastq.gz
    """
}

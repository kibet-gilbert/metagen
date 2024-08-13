process METASPADES {
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
        -1 ${nohost_reads[0]} \\
        -2 ${nohost_reads[1]} \\
        -o .
    mv assembly_graph_with_scaffolds.gfa SPAdes-${sample_id}_graph.gfa         
    mv scaffolds.fasta SPAdes-${sample_id}_scaffolds.fasta                     
    mv contigs.fasta SPAdes-${sample_id}_contigs.fasta                         
    mv spades.log SPAdes-${sample_id}.log                                      
    """
}


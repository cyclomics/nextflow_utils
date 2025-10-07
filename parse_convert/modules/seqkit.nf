process MergeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path fasta1
        path fasta2
    
    output:
        path "${fasta1.simpleName}_${fasta2.simpleName}.fasta"
    
    script:
        """
        cat $fasta1 > ${fasta1.simpleName}_${fasta2.simpleName}.fasta
        cat $fasta2 >> ${fasta1.simpleName}_${fasta2.simpleName}.fasta
        """

}

process SplitReadFiles{
    label 'many_cpu_medium'
    
    input:
        tuple val(sample_id), val(fq_id), path(fastq)

    output:
        tuple val(sample_id), val(fq_id), path ("split/${fastq.simpleName}_part_*.fastq*")

    script:
        // Split the fastq file into smaller files if they are more than reads_per_fq read
        // --by-size-prefix prevents . in the filenames, which can cause issues like filename collisions with .simpleName
        """
        seqkit split --by-size ${params.splitting.reads_per_fq} --by-size-prefix ${fastq.simpleName}_part_ -O split $fastq
        """
}


process Extract5PrimeFasta {
    input:
        path fasta
        val length

    output:
        path('*.fastq') 

    script:
        """
        seqtk trimfq -L $length $fasta 
        """
}


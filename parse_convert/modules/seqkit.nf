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


process FilterShortReads{
    label 'many_cpu_medium'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    
    input:
        path(fastq)

    output:
        path ("${fastq.simpleName}_filtered.fastq")

    script:
        """
        seqkit seq -m ${params.filtering.minimun_raw_length} $fastq > "${fastq.simpleName}_filtered.fastq"
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


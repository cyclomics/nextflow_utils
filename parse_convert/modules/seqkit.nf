
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


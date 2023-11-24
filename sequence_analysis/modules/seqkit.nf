
process Protocol{
    // Example of IO for each process
    input:
        tuple val(sample), val(ID), path(fq)
    output:
        tuple val(sample), val(ID), path("*.fq")
    script:
        """
        cat $fq > ${sample}_${ID}_protocol.fq
        """
}

process FilterShortReads{   
    input:
        tuple val(sample), val(ID), path(fq)

    output:
        tuple val(sample), val(ID), path(fq)

    script:
        """
        seqkit seq -m ${params.filtering.minimun_raw_length} $fq > "${fq.simpleName}_filtered.fastq"
        """
}
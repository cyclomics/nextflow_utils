
process FilterShortReads{   
    input:
        tuple val(sample), val(ID), path(fq)

    output:
        tuple val(sample), val(ID), path("${fq.simpleName}_filtered.fastq")

    script:
        """
        seqkit seq -m ${params.filtering.minimun_raw_length} $fq > "${fq.simpleName}_filtered.fastq"
        """
}
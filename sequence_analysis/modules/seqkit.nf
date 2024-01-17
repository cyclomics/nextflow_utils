
process FilterShortReads{
    // Remove all reads shorter than filtering.minimun_raw_length
    input:
        tuple val(sample), val(ID), path(fq)

    output:
        tuple val(sample), val(ID), path("${fq.simpleName}_filtered.fastq")

    script:
        """
        seqkit seq -m ${params.filtering.minimun_raw_length} $fq > "${fq.simpleName}_filtered.fastq"
        """
}


process FilterShortReads{
    // Remove all reads shorter than filtering.minimun_raw_length
    input:
        path(fq)

    output:
        path("${fq.simpleName}_filtered.fastq")

    script:
        """
        seqkit seq -m ${params.filtering.minimun_raw_length} $fq > "${fq.simpleName}_filtered.fastq"
        """
}

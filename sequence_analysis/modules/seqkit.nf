
process FilterShortReads{
    // Remove all reads shorter than filtering.minimun_raw_length
    
    input:
        tuple val(sample_id), val(file_id), path(fq)

    output:
        tuple val(sample_id), val("${fq.simpleName}_filtered"), path("${fq.simpleName}_filtered.fastq")

    script:
        """
        seqkit seq -m ${params.filtering.minimun_raw_length} $fq > "${fq.simpleName}_filtered.fastq"
        """
}

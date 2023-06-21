
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

process TrimPolyATails{
    input:
        tuple val(sample), val(ID), path(fq)
    output:
        tuple val(sample), val(ID), path("*polyatrim.fq")
    script:
        """
        cutadapt -ploy-a -o ${fq.simpleName}.polyatrim.fastq.gz $fq
        """
}


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

process Cut5P3PAdapter {
    input:
        tuple val(sample), val(ID), path(fq)
        val(adapter_5p)
        val(adapter_3p)
    output:
        tuple val(ID), path("${fq.simpleName}_cutadapt.fq.gz")
    script:
    """
    cutadapt -a $adapter_3p -g $adapter_5p -o ${fq.simpleName}_cutadapt.fq.gz $fq
    """
}

process Extract5PrimeFasta{
    input:
        path(fa)
        val(requested_length)
    output:
        stdout
    script:
    // see: https://bioinf.shenwei.me/seqkit/usage/#subseq
    // Dont forget to trim the newline char
    """
    seqkit head -n 1 $fa | seqkit subseq -r 1:$requested_length | seqkit seq --seq | tr -d '\n'
    """
}
process Extract3PrimeRemainderFasta{
    input:
        path(fa)
        val(requested_length)
    output:
        stdout
    script:
    // see: https://bioinf.shenwei.me/seqkit/usage/#subseq
    """
    # seqkit stats $fa
    READ_COUNT=\$(seqkit stats $fa | awk ' NR==1 {next} {print \$4}')
    READ_LENGTH=\$(seqkit stats $fa | awk ' NR==1 {next} {print \$6}')
    # echo \$READ_COUNT
    # echo \$READ_LENGTH
    READ_CUT=\$(expr \$READ_LENGTH - $requested_length)
    # echo \$READ_CUT
    seqkit head -n 1 $fa | seqkit subseq -r -\$READ_CUT:-1 | seqkit seq --seq | tr -d '\n'

    """
}

process ExtractFullFasta{
    input:
        path(fa)
    output:
        stdout
    script:
    // see: https://bioinf.shenwei.me/seqkit/usage/#subseq
    // Dont forget to trim the newline char
    """
    seqkit head -n 1 $fa | seqkit seq --seq | tr -d '\n'
    """
}

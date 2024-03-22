
process Protocol{
    // Example of IO for each process in this modules folder
    input:
        tuple val(sample), val(ID), path(fq)
    output:
        tuple val(sample), val(ID), path("*.fq")
    script:
        """
        cat $fq > ${sample}_${ID}_protocol.fq
        """
}


process SummerizeReadsStdout{
    // Write stats to stdout, does not alter the input fq in any way.
    input:
        tuple val(sample_ID), val(file_ID), file(samples)

    output:
        stdout

    script:
        """
        echo Summarizing reads for__: $sample_ID
        echo Example file name______: ${samples.first()}
        echo number of files________: ${samples.size()}
        echo combined statistics for the sample:
        seqkit seq \$(find . -name "*.fq*" -o -name "*.fastq*") | seqkit stats -ab
        """
}


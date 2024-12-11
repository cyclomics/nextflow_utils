process SplitReadsOnAdapterSequence {
    // https://github.com/nanoporetech/duplex-tools/blob/master/fillet.md
    label 'many_low_cpu_low_mem'

    input:
        tuple val(sample_id), val(file_id), path(fq)

    output:
        tuple val(sample_id), val(file_id), path("results/${fq.simpleName}_split.fastq.gz")
        
    script:
        """
        duplex_tools split_on_adapter . results/ Native 
        """
} 
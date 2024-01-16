#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SamtoolsMergeBams{
    //  merge n number of bams into one
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    label 'many_cpu_medium'

    input:
        val(X)
        file(bam_in)

    output:
        tuple val(X), path("${X}.merged.bam"), path("${X}.merged.bam.bai")
    
    script:
    """
    ls
    samtools merge -p -c -O bam ${X}.merged.bam \$(find . -name '*.bam')
    samtools index ${X}.merged.bam
    """
}
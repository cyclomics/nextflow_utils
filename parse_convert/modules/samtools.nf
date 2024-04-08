#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SamtoolsMergeBams{
    //  merge n number of bams into one
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    label 'many_cpu_medium'

    input:
        val(sample_id)
        path(bam_in)

    output:
        tuple val(sample_id), path("${sample_id}.merged.bam")
    
    script:
    """
    samtools merge -p -c -O bam ${sample_id}.merged.bam \$(find . -name '*.bam')
    samtools index ${sample_id}.merged.bam
    """
}
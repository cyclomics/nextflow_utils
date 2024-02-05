#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SamtoolsMergeBams{
    //  merge n number of bams into one
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    label 'many_cpu_medium'

    input:
        tuple val(sample_id), val(file_id), file(bam_in)

    output:
        tuple val(sample_id), val("${bam_in.simpleName}.merged"), path("${bam_in.simpleName}.merged.bam"), path("${bam_in.simpleName}.merged.bam.bai")
    
    script:
    """
    ls
    samtools merge -p -c -O bam ${bam_in.simpleName}.merged.bam \$(find . -name '*.bam')
    samtools index ${bam_in.simpleName}.merged.bam
    """
}
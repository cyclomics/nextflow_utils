#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process Cygnus{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/consensus", mode: 'copy'
    label 'many_low_cpu_tiny_mem'

    input:
        tuple val(X), path(input_folder)

    output:
        tuple val(X), path("${bam.SimpleName}.consensus.fastq")

    script:
        
        """
        python $params.cygnus_location $input_folder $output_folder --threads=$threads
        """
    }
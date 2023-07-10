#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process Cycas {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/consensus", mode: 'copy'
    label 'many_low_cpu_tiny_mem'

    input:
        tuple val(X), path(bam), path(bai)

    output:
        tuple val(X), path("${bam.SimpleName}.consensus.fastq")

    script:
        
        """
        mkdir plots
        python $params.cycas_location consensus --bam-file $bam --output ${bam.SimpleName}.consensus.fastq --metadata-json  ${bam.SimpleName}.metadata.json
        """
    }
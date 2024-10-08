#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process Cycas {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // publishDir "${params.output_dir}/consensus", mode: 'copy'
    label 'many_low_cpu_tiny_mem'
    container = "cyclomics/cycas:0.5.2"
    cpus 1
    memory = '2G'

    input:
        tuple val(sample_id), val(file_id), path(bam), path(bai)

    output:
        tuple val(sample_id), val(file_id), path("${bam.simpleName}.consensus.fastq"), path("${bam.simpleName}.metadata.json")

    script:
        """
        mkdir -p plots
        python $params.cycas_location consensus --bam-file $bam --output ${bam.simpleName}.consensus.fastq --metadata-json  ${bam.simpleName}.metadata.json
        """
    }

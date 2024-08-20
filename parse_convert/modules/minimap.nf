#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process Minimap2Index{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'few_memory_intensive'
    
     input:
        path(reference_genome)
    
    output:
        path("${reference_genome.simpleName}.mmi")

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} -d ${reference_genome.simpleName}.mmi $reference_genome
        """
}

process Minimap2AlignAdaptive{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'minimap_large'

    //  apply at least 1 Gb of memory to the process, otherwise we apply two-four times the size of the reference genome mmi
    memory = {reference_genome.size() > 1_000_000_000 ? Math.round(reference_genome.size()*1.8 * task.attempt) : "4GB"* task.attempt}
    // memory "32 GB"
    // small jobs get 4 cores, big ones 8
    cpus (params.economy_mode == true ? 2 :{reference_genome.size() < 500_000_000 ? 4 : 8 })
    
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id), val(file_id), path(fastq)
        path(reference_genome)
    
    output:
        tuple val(sample_id), val("${fastq.simpleName}"), path("${fastq.simpleName}.bam")

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} $reference_genome $fastq > tmp.sam 
        samtools sort -o ${fastq.simpleName}.bam tmp.sam
        rm tmp.sam
        """
}

process Minimap2AlignAdaptiveParameterized{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    
    // *1.8 for T2T gives 13.5 for the first try and 27 for the second
    // Thus fitting within 32GB systems (which have less than 32 in reality)
    // Small refs below 1GB get 4GB per task retry (4,8,12)
    memory = {reference_genome.size() > 1_000_000_000 ? Math.round(reference_genome.size()*1.8 * task.attempt) : "4GB"* task.attempt}
    cpus params.economy_mode == true ? 2 : 7

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        // each path(fastq)
        tuple val(sample_id), val(file_id), path(fastq)
        path(reference_genome)


    output:
        tuple val(sample_id), val("${fastq.simpleName}"), path("${fastq.simpleName}.bam") 

    script:
    // Lower parameters to increase data available to cycas
        """
        minimap2 -ax map-ont -t ${task.cpus} -m ${params.minimap2parameterized.min_chain_score} -n ${params.minimap2parameterized.min_chain_count} -s ${params.minimap2parameterized.min_peak_aln_score} $reference_genome $fastq | samtools sort -o ${fastq.simpleName}.bam
        """
}

process Minimap2Align{    
    // Use standard Minimap2 parameters for alignment, also works with .mmi files.
    cpus params.economy_mode == true ? 2 : 7
    
    // *1.8 for T2T gives 13.5 for the first try and 27 for the second
    // Thus fitting within 32GB systems (which have less than 32 in reality)
    // Small refs below 1GB get 4GB per task retry (4,8,12)
    memory = {reference_genome.size() > 1_000_000_000 ? Math.round(reference_genome.size()*1.8 * task.attempt) : "4GB"* task.attempt}
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id), val(fq_id), path(fq)
        path(reference_genome)
    
    output:
        tuple val(sample_id), val(fq_id), path("${fq_id}.bam")

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} $reference_genome $fq > tmp.sam 
        samtools sort -o ${fq_id}.bam tmp.sam
        rm tmp.sam
        """
}

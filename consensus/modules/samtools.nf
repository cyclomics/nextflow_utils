#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process SamtoolsIndexWithID{
    label 'many_cpu_medium'

    input:
        tuple val(x), path(input_bam) 

    output:
        tuple val(x), path(input_bam), path("*.bai")

    script:
        """
        samtools index $input_bam
        """
}

process PrimaryMappedFilter{
    // Sort, convert and index 
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam_in), path(bai_in)

    output:
        tuple val(X), path("${X}.primary_mapped.bam"), path("${X}.primary_mapped.bam.bai")

    script:
        """
        samtools view -b -F 256 $bam_in > ${X}.primary_mapped.bam
        samtools index ${X}.primary_mapped.bam
        """
}

process MapqAndNMFilter{
    // Sort, convert and index 
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam_in), path(bai_in)

    output:
        tuple val(X), path("${X}.NM_50_mapq_20.bam"), path("${X}.NM_50_mapq_20.bam.bai")

    script:
        """
        samtools view -b -o ${X}.NM_50_mapq_20.bam $bam_in --input-fmt-option 'filter=[NM]<50 && mapq >20'
        samtools index ${X}.NM_50_mapq_20.bam
        """
}
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process Tidehunter53QualTable {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'quay.io/biocontainers/tidehunter:1.5.3--h2e03b76_0'
    label 'many_cpu_medium'

    input:
        tuple path(fasta), path(prime_3_fasta), path(prime_5_fasta)

    output:
        tuple val("${fasta.baseName}"), path("*consensus.tsv")

    script:
        """
        echo "$params.tidehunter.headerlinesQual" > ${fasta.SimpleName}.consensus.tsv
        TideHunter \
        --out-fmt 4  \
        --kmer-length 16 \
        --longest  \
        --thread ${task.cpus} \
        --five-prime $prime_5_fasta \
        --three-prime $prime_3_fasta \
        --min-period $params.tidehunter.minimum_period \
        --min-len $params.tidehunter.minimum_length \
        --min-copy $params.tidehunter.minimum_copy \
        -a $params.tidehunter.minimum_match_ratio \
        $fasta >> ${fasta.SimpleName}.consensus.tsv
        """
}

process TideHunterFilterTableStartpos {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(tidehuntertable)
    
    output:
        tuple val(X), path("${tidehuntertable.baseName}.startposfilter.tsv")

    script:
        """
        awk 'BEGIN {getline; print}{if (\$5 < 100) {print}}' $tidehuntertable > ${tidehuntertable.baseName}.startposfilter.tsv
        """
}

process TideHunterQualTableToFastq {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(tidehuntertable)

    
    output:
        tuple val(X), path("${tidehuntertable.SimpleName}.fastq")

    script:
        //  skip first line (header), go on to convert all lines to fastq entries
        """
        awk 'BEGIN {getline}{print "@"\$1"\\n"\$11"\\n+\\n"\$12}' $tidehuntertable > ${tidehuntertable.SimpleName}.fastq
        """
}
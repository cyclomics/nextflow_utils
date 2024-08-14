#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process Cygnus {
    container = "cyclomics/cygnus:0.6.5"
    cpus 1
    memory = '2G'

    input:
        tuple val(sample_id), val(fq_id), path(fq)
    output:
        tuple val(sample_id), val(fq_id), path("*_cygnus.fq.gz")
    script:
        // Check if python file for cygnus exist in expected location, execute it when found 
        """
        if [[ -f "${params.cygnus_location}" ]]; then
            echo "Executing Cygnus with input file: $fq"
            python $params.cygnus_location $fq ${fq.SimpleName}_cygnus.fq.gz
        else
            echo "Cygnus submodule not found at ${params.cygnus_location}, please make sure the software is available."
            exit 1
        fi
        """
}

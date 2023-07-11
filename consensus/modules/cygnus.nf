#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process Cygnus {
    input:
        tuple val(sample_id), val(fq_id), path(fq)
    output:
        tuple val(sample_id), val(fq_id), path("*.cygnus.fq.gz")
    script:
        // Check if python file for cygnus exist in expected location, execute it when found 
        """
        if [[ -f "${params.cygnus_location}" ]]; then
            touch ${fq.SimpleName}.cygnus.rejects.fq.gz
            python $params.cygnus_location $fq ${fq.SimpleName}.cygnus.fq.gz --rejects_file ${fq.SimpleName}.cygnus.rejects.fq.gz
        else
            echo "Cygnus submodule not found, please add it using git submodules to this pipeline."
            exit 1
        fi
        """
}

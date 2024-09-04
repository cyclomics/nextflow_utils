
process Cygnus2 {
    container = "cyclomics/cygnus2:0.0.1"
    cpus 1
    memory = '2G'

    input:
        tuple val(sample_id), val(fq_id), path(fq)
    output:
        tuple val(sample_id), val(fq_id), path("*_cygnus.fq.gz")
    script:
        // Check if python file for cygnus exist in expected location, execute it when found 
        """
        if [[ -f "${params.cygnus2_location}" ]]; then
            echo "Executing Cygnus with input file: $fq"
            python $params.cygnus2_location run $fq --output_path ${fq.SimpleName}_cygnus.fq.gz
        else
            echo "Cygnus2 cmd not found at ${params.cygnus2_location}, please make sure the software is available."
            exit 1
        fi
        """
}

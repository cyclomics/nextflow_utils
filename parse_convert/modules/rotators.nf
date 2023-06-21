process RotateBySequence{
     input:
        tuple val(sample), val(ID), path(fq)
    output:
        tuple val(sample), val(ID), path("*.rotated.fq")
    script:
    """
    if [[ -f "${params.rotator_location}/rotate_by_sequence.py" ]]; then
        python ${params.rotator_location}/rotate_by_sequence.py
    else
        echo "Rotator submodule not found, please add it using git submodules to this pipeline."
        exit 1
    fi
    """
}
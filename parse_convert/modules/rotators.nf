process RotateBySequence {
    input:
        tuple val(ID), path(fq)
    output:
        tuple val(ID), path("*.rotated.fq")
    script:
        """
        if [[ -f "${params.rotator_location}/rotate_by_sequence.py" ]]; then
            python ${params.rotator_location}/rotate_by_sequence.py $fq ${fq.SimpleName}.rotated.fq
        else
            echo "Rotator submodule not found, please add it using git submodules to this pipeline."
            exit 1
        fi
    """
}
process RotateByAlignment {
    container = "cyclomics/rotators:0.1.1"
    cpus 1
    memory "2GB"
    input:
        tuple val(sample_id), val(fq_id), path(bam)
    output:
        tuple val(sample_id), val(fq_id), path("${bam.SimpleName}_rotated.fq")
    script:
        """
        if [[ -f "${params.rotator_location}/rotate_bam_by_alignment.py" ]]; then
            python ${params.rotator_location}/rotate_bam_by_alignment.py -i $bam -o ${bam.SimpleName}_rotated.fq
        else
            echo "Rotator submodule not found, please add it using git submodules to this pipeline."
            exit 1
        fi
        """
}
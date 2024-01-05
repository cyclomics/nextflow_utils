nextflow_process {

    name "Test Process SamtoolsIndexWithID"
    script "consensus/modules/samtools.nf"
    process "SamtoolsIndexWithID"
    
    
    test("Samtools BAM indexing") {

        when {
            process {
                """
                sample_id = "test_sample"
                bam_file = file("$baseDir/example_datasets/example_native_reads/bam/FAW79986_pass_dbcafcd9_40623079_0_filtered.bam")

                input[0] = Channel.of([sample_id, bam_file])
                """
            }
        }
        then {
            assert process.success
            
            def bai_filename = "FAW79986_pass_dbcafcd9_40623079_0_filtered.bam.bai"
            with(process.out[0]) {
                assert path(get(0)[2]).getFileName().toString() == bai_filename
                assert path(get(0)[2]).size() > 0
            }
        }
    }
}

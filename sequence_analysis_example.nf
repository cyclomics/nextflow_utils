

include {
    TrimPolyATails
} from './sequence_analysis/modules/cutadapt'

// Cyclomics Nextflow utils primary analysis example

params.input_read_regex = "example_datasets/example_consensus_reads/reads/*.fastq"
params.sample_id = " undefined"

workflow {
    read_fastq = Channel.fromPath(params.input_read_regex, checkIfExists: true) | 
        map(x -> [params.sample_id, x.Name, x])

    read_fastq.view()

    adapter_trimmed_fastq = TrimPolyATails(read_fastq)

    adapter_trimmed_fastq.view()
}
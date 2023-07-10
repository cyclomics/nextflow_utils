nextflow.enable.dsl = 2

include {
    FilterShortReads
} from "./modules/seqkit"

include {
    AnnotateBamXTags
} from "./modules/annotate_bam.nf"

process Protocol{
    // Example of IO for each process
    input:
        tuple val(sample), val(ID), path(fq)
    output:
        tuple val(sample), val(ID), path("*.fq")
    script:
        """
        cat $fq > ${sample}_${ID}_protocol.fq
        """
}

workflow FilterWithAdapterDetection{
    take:
        read_fq_ch

    main:
        if (params.split_on_adapter == true) {
            fastq = SplitReadsOnAdapterSequence(read_fq_ch)
        }
        else {
            fastq = read_fq_ch
        }
        FilterShortReads(fastq)
        
    emit:
        FilterShortReads.out
}

workflow AnnotateBam{
    take:
        reads
        sequencing_summary

    main:
        AnnotateBamXTags(reads, sequencing_summary)
    emit:
        AnnotateBamXTags.out
}

workflow FilterBam{
    take:
        annotated_bam
        minimun_repeat_count

    main:
        BamTagFilter(annotated_bam, 'YM', minimun_repeat_count)
    emit:
        BamTagFilter.out
}
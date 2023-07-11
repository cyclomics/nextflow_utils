#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
    FilterShortReads
} from "./modules/seqkit"

include {
    AnnotateBamXTags
} from "./modules/annotate_bam"

include {
    BamTagFilter
} from "./modules/samtools"

workflow FilterWithAdapterDetection {
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

workflow AnnotateBam {
    take:
        bam

    main:
        seqsum = Channel.fromPath(params.sequencing_summary_path, checkIfExists: true)
        AnnotateBamXTags(bam, seqsum)
        
    emit:
        AnnotateBamXTags.out
}

workflow FilterBam {
    take:
        bam

    main:
        BamTagFilter(bam, 'YM', params.min_repeat_count)

    emit:
        BamTagFilter.out
}
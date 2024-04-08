#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
    MergeFasta
} from "./modules/seqkit"

include {
    Minimap2AlignAdaptive
    Minimap2Index as IndexCombined
} from "./modules/minimap.nf"

include {
    SamtoolsMergeBams
    // BamTagFilter
} from "./modules/samtools.nf"

// include {
//     BwaIndex
//     BwaMemSorted
//     // SamtoolsMergeBams
// } from "./modules/bwa.nf"

workflow PrepareGenome {
    take:
        reference
        reference_name
        backbones_fasta

    main:
        if (reference_name.endsWith('.txt')) {
            println("txt reference, not implemented. Exiting...")
            genome = "missing"
            exit(1)
        }
        else if (reference_name.endsWith('.gz')) {
            println("gzipped reference, not implemented. Exiting...")
            genome = "missing"
            exit(1)
        }
        else {
            genome = reference
        }

        MergeFasta(genome, backbones_fasta)
        IndexCombined(MergeFasta.out)
        
    emit:
        mmi_combi = IndexCombined.out
        fasta_combi = MergeFasta.out
}

workflow Minimap2Align {
    take:
        reads
        reference

    main:
        id = reads.first().map( it -> it[0])
        id = id.map(it -> it.split('_')[0])
        Minimap2AlignAdaptive(reads, reference)
        // id = Minimap2AlignAdaptive.out.map(it -> it[0].split('_')[0])
        bams = Minimap2AlignAdaptive.out.map(it -> it[2]).collect()
        SamtoolsMergeBams(id, bams)

    emit:
        bam = SamtoolsMergeBams.out
}

// workflow BWAAlign{
//     take:
//         reads
//         reference

//     main:
//         bwa_idx = file("${params.reference}.{,amb,ann,bwt,pac,sa}")
//         bwa_index_file_count = 5
//         // We do a smaller than since there might be a .fai file as well!
//         if (bwa_idx.size < bwa_index_file_count) {
//             println "==================================="
//             println "Warning! BWA index files are missing for the reference genome, This will slowdown execution in a major way."
//             println "==================================="
//             println ""
//             println ""
//             bwa_idx = BwaIndex(reference)
//         }
        
//         id = reads.first().map(it -> it[0])
//         id = id.map(it -> it.split('_')[0])
//         reads_fastq = reads.map(it -> it[1])
//         BwaMemSorted(reads_fastq, reference, bwa_idx.collect())
//         SamtoolsMergeBams(id, BwaMemSorted.out.collect())

//     emit:
//         bam = SamtoolsMergeBams.out
// }
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
    Cygnus
} from "./modules/cygnus"

include {
    Cycas
} from "./modules/cycas"

include {
    Extract5PrimeFasta
    Extract3PrimeFasta
    ExtractSpecificRead
} from "./modules/seqkit"

include{
    Tidehunter53QualTable
    TideHunterFilterTableStartpos
    TideHunterQualTableToFastq
} from "./modules/tidehunter"

include {
    RotateByAlignment
    RotateBySequence
} from "../parse_convert/modules/rotators"

include {
    Minimap2AlignAdaptiveParameterized
    Minimap2Align
    Minimap2Index
} from "../parse_convert/modules/minimap"

include {
    SamtoolsIndexWithID
    PrimaryMappedFilter
    MapqAndNMFilter
} from "./modules/samtools"


workflow CygnusConsensus {
    take:
        reads_fastq

    main:
        Cygnus(reads)
        RotateBySequence(Cygnus.out)

    emit:
        fastq = RotateBySequence.out

}

workflow CygnusAlignedConsensus {
    take:
        reads_fastq
        reference

    main:
        // 1. Create reads
        Cygnus(reads_fastq)
        consensus = Cygnus.out
        // 2. Align the reads
        Minimap2Align(consensus, reference)
        // 3. Rotate by alignment
        RotateByAlignment(Minimap2Align.out)

    emit:
        fastq = RotateByAlignment.out

}
workflow CygnusPrimedConsensus {
    take:
        reads_fastq

    main:
        Cygnus(reads)

    emit:
        fastq = Cygnus.out

}

workflow CycasConsensus {
    take:
        read_fastq
        reference_genome

    main:
        // Its the callers responsibility to make sure reference_genome is a value channel
        Minimap2AlignAdaptiveParameterized(read_fastq, reference_genome)
        SamtoolsIndexWithID(Minimap2AlignAdaptiveParameterized.out)
        PrimaryMappedFilter(SamtoolsIndexWithID.out)
        MapqAndNMFilter(PrimaryMappedFilter.out)
        Cycas(MapqAndNMFilter.out)

    emit:
        Cycas.out
        // fastq = Cycas.out.map(it -> it.take(3))
        // json = Cycas.out.map(it -> tuple(it[0], it[3]))
}

workflow TidehunterConsensus {
    take:
        reads_fastq
        reference_genome
        backbone_fasta

    main: 
        ExtractSpecificRead(backbone_fasta, params.backbone_name)
        Extract5PrimeFasta(ExtractSpecificRead.out, params.tidehunter.primer_length)
        Extract3PrimeFasta(ExtractSpecificRead.out, params.tidehunter.primer_length)

        Tidehunter53QualTable(read_fastq.combine(Extract5PrimeFasta.out).combine(Extract3PrimeFasta.out))
        TideHunterFilterTableStartpos(Tidehunter53QualTable.out)
        TideHunterQualTableToFastq(TideHunterFilterTableStartpos.out)
    
    emit:
        fastq = TideHunterQualTableToFastq.out
        // json = TideHunterQualJsonMerge.out
}
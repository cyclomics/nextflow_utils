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
    RotateBySequence
} from "../parse_convert/modules/rotators"

include {
    Minimap2AlignAdaptiveParameterized
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

workflow CycasConsensus {
    take:
        reads_fastq
        reference_genome

    main:
        Minimap2AlignAdaptiveParameterized(reads_fastq, reference_genome)
        SamtoolsIndexWithID(Minimap2AlignAdaptiveParameterized.out)
        PrimaryMappedFilter(SamtoolsIndexWithID.out)
        MapqAndNMFilter(PrimaryMappedFilter.out)
        Cycas(MapqAndNMFilter.out)

    emit:
        fastq = Cycas.out.map(it -> it.take(2))
        json = Cycas.out.map(it -> tuple(it[0], it[2]))
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
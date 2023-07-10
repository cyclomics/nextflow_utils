nextflow.enable.dsl = 2

include {
    Cygnus
} from "./modules/cygnus"

include {
    Cycas
} from "./modules/Cycas"


include {
    RotateBySequence
} from "../parse_convert/modules/rotators"

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
        backbone_fasta
    main:
        Minimap2AlignAdaptiveParameterized(reads_fastq, reference_genome)
        SamtoolsIndexWithID(Minimap2AlignAdaptiveParameterized.out)
        PrimaryMappedFilter(SamtoolsIndexWithID.out)
        MapqAndNMFilter(PrimaryMappedFilter.out)
        Cycas(MapqAndNMFilter.out)

    emit:
        fastq = Cycas.out
        // take the id and json
        // json = Cycas.out.map( it -> tuple(it[0], it[2]))
        // split_bam = Minimap2AlignAdaptiveParameterized.out
        // split_bam_filtered = MapqAndNMFilter.out
}
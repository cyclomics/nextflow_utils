nextflow.enable.dsl = 2

include {
    PrimaryMappedFilter
    SamtoolsIndexWithID
    MapqAndNMFilter
} from "./modules/samtools"

include {
    Extract5PrimeFasta
    Extract3PrimeFasta
    ExtractSpecificRead
} from "./modules/seqkit"

include {
    Tidehunter53QualTable
    TideHunterFilterTableStartpos
    TideHunterQualTableToFastq
    TideHunterQualTableToJson
    TideHunterQualJsonMerge
} from "./modules/tidehunter"

include {
    Minimap2AlignAdaptiveParameterized
    Minimap2AlignAdaptiveParameterized as MinimapForSplitMap
} from "./modules/minimap"

include {
    Cycas
} from "./modules/cycas"

// workflow CygnusConsensus {}

workflow CycasConsensus {
    take:
        read_fastq
        reference_genome
        backbone_fasta
    main:
        Minimap2AlignAdaptiveParameterized(read_fastq, reference_genome)
        SamtoolsIndexWithID(Minimap2AlignAdaptiveParameterized.out)
        PrimaryMappedFilter(SamtoolsIndexWithID.out)
        MapqAndNMFilter(PrimaryMappedFilter.out)
        Cycas(MapqAndNMFilter.out)

    emit:
        // take the first 2 elements aka id and fastq
        fastq = Cycas.out.map(it -> it.take(2))
        // take the id and json
        json = Cycas.out.map(it -> tuple(it[0], it[2]))
        split_bam = Minimap2AlignAdaptiveParameterized.out
        split_bam_filtered = MapqAndNMFilter.out
}

workflow TidehunterBackBoneQual {
    // TidehunterBackBoneQual takes the backbones and runs tidehunter whilst providing the 3 and 5 prime regions to tidehunter,
    // TODO: rotate?

    take:
        read_fastq
        reference_genome
        backbone_fasta
        backbone_primer_len
        backbone_name
        reference_genome_mmi

    main:
        // get backbones
        ExtractSpecificRead(backbone_fasta, backbone_name)
        Extract5PrimeFasta(ExtractSpecificRead.out, backbone_primer_len)
        Extract3PrimeFasta(ExtractSpecificRead.out, backbone_primer_len)

        Tidehunter53QualTable(read_fastq.combine(Extract5PrimeFasta.out).combine(Extract3PrimeFasta.out))
        TideHunterFilterTableStartpos(Tidehunter53QualTable.out)
        TideHunterQualTableToFastq(TideHunterFilterTableStartpos.out)
        TideHunterQualTableToJson(TideHunterFilterTableStartpos.out)

        id = TideHunterQualTableToJson.out.first()map(it -> it[0])
        id = id.map(it -> it.split('_')[0])
        jsons = TideHunterQualTableToJson.out.map(it -> it[1]).collect()
        TideHunterQualJsonMerge(id, jsons)
        MinimapForSplitMap(read_fastq, reference_genome_mmi)

    emit:
        fastq = TideHunterQualTableToFastq.out
        json = TideHunterQualJsonMerge.out
        split_bam = MinimapForSplitMap.out
        split_bam_filtered = MinimapForSplitMap.out
}
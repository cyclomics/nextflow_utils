nextflow.enable.dsl = 2

include {
    CountFastqInfo as FastqInfoRaw
    CountFastqInfo as FastqInfoConsensus
} from "./modules/seqkit"

include {
    SamtoolsQuickcheck
    SamtoolsFlagstats
    SamtoolsIdxStats
    SamtoolsMergeBams as SamtoolsMergeBams
    SamtoolsMergeBams as SamtoolsMergeBamsFiltered
    FindRegionOfInterest
} from "./modules/samtools"

include {
    PerbaseBaseDepth as PerbaseBaseDepthSplit
    PerbaseBaseDepth as PerbaseBaseDepthConsensus
} from "./modules/perbase"

include {
    PlotReadStructure
    PlotFastqsQUalAndLength as PlotRawFastqHist
    PlotFastqsQUalAndLength as PlotFilteredHist
    PlotFastqsQUalAndLength as PlotConFastqHist
    PlotQScores
    PlotMetadataStats
    PlotReport
} from "./modules/bin"

workflow PostQC {
    take:
        reference_fasta
        fastq_raw
        fastq_filtered
        split_bam
        split_bam_filtered
        fastq_consensus
        read_info
        consensus_bam

    main:
        // dont add the ID to the process
        first_fq = fastq_raw.first()
        id = first_fq.simpleName
        extension = first_fq.getExtension()
        FastqInfoRaw(fastq_raw.collect(),'raw')
        PlotRawFastqHist(fastq_raw.collect(), extension, id + "raw", '"raw fastq info"')
        
        first_fq = fastq_filtered.first()
        id = first_fq.simpleName
        extension = first_fq.getExtension()
        PlotFilteredHist(fastq_filtered.collect(),extension, id + "filtered", '"filtered fastq info"')

        first_fq = fastq_consensus.first()
        id = first_fq.map(it -> it[0])
        extension = first_fq.map(it -> it[1]).getExtension()

        FastqInfoConsensus(fastq_consensus.map(it -> it[1]).collect(), 'consensus')
        PlotConFastqHist(fastq_consensus.map(it -> it[1]).collect(),extension, id + "consensus", '"consensus fastq info"')

        merged_split_bam = SamtoolsMergeBams('splibams_merged', split_bam.collect())
        PlotReadStructure(merged_split_bam)
        SamtoolsQuickcheck(consensus_bam)
        SamtoolsIdxStats(consensus_bam)
        meta_data = read_info.map(it -> it[1]).collect()
        PlotMetadataStats(meta_data)

        roi = FindRegionOfInterest(consensus_bam)
      
        PerbaseBaseDepthSplit(merged_split_bam.combine(reference_fasta), roi, 'split.tsv')
        PerbaseBaseDepthConsensus(consensus_bam.combine(reference_fasta), roi, 'consensus.tsv')

        PlotQScores(PerbaseBaseDepthSplit.out, PerbaseBaseDepthConsensus.out)

        merged_split_bam_filtered = SamtoolsMergeBamsFiltered('splibams_filtered_merged',split_bam_filtered.collect())
        SamtoolsFlagstats(merged_split_bam_filtered)
        
        if (params.report == true) {
            PlotReport(
                PlotRawFastqHist.out.combine(
                PlotFilteredHist.out).combine(
                PlotConFastqHist.out).combine(
                PlotReadStructure.out).combine(
                PlotQScores.out).combine(
                PlotMetadataStats.out).combine(
                PasteVariantTable.out).combine(
                SamtoolsFlagstats.out).combine(
                SamtoolsIdxStats.out
                )
            )
        }
}
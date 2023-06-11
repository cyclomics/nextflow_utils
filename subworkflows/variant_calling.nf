nextflow.enable.dsl=2

include {
    Freebayes
} from "./modules/freebayes"

include {
    VcfToBed
} from "./modules/bedops"

include {
    FindVariants
    FilterVariants
    MergeNoisyVCF
    MergeFilteredVCF
    AnnotateVCF
} from "./modules/variant_validate"

include {
    FindRegionOfInterest
} from "./modules/samtools"

include {
    PerbaseBaseDepth as PerbaseBaseDepthConsensus
} from "./modules/perbase"


workflow ProcessTargetRegions{
    take:
        variant_file_name
        reads_aligned
    
    main:
        variant_file = Channel.fromPath(params.region_file, checkIfExists: false)

        if (variant_file_name == 'auto') {
            positions = FindRegionOfInterest(reads_aligned)
        }
        else if( variant_file_name.endsWith('.bed') ) {
            positions = variant_file
        }
        else if( variant_file_name.endsWith('.vcf') ) {
            positions = VcfToBed(variant_file)
        }
        else {
            positions = variant_file.map(it -> tuple(it.SimpleName, it))
        }
    
    emit:
        positions

}

workflow FreebayesSimple{
    take:
        reads_aligned
        positions
        reference

    main:
        Freebayes(reads_aligned.combine(reference), positions)
        AnnotateVCF(Freebayes.out)

    emit:
        locations = Freebayes.out
        variants = AnnotateVCF.out
}


workflow ValidatePosibleVariantLocations{
    // Allow to determine VAF for given genomic positions in both bed and vcf format
    take:
        reads_aligned
        positions
        reference

    main:
        FindVariants(reference, reads_aligned, positions)
        PerbaseBaseDepthConsensus(reads_aligned.combine(reference), positions, 'consensus.tsv')
        FilterVariants(FindVariants.out.combine(PerbaseBaseDepthConsensus.out))
        MergeNoisyVCF(FindVariants.out)
        MergeFilteredVCF(FilterVariants.out)
        AnnotateVCF(MergeFilteredVCF.out)

    emit:
        locations = MergeNoisyVCF.out
        variants = AnnotateVCF.out
}
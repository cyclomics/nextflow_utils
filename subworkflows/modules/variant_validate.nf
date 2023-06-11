nextflow.enable.dsl=2

process FindVariants{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'max_performance'

    input:
        path(reference_genome)
        tuple val(X), path(bam), path(bai)
        tuple val(X), path(validation_bed)

    output:
        tuple path("${bam.simpleName}.snp.vcf"), path("${bam.simpleName}.indel.vcf")
    
    script:
        // We sleep and access the reference genome, since in some rare cases the file needs accessing to 
        // not cause issues in the python code. 
        """
        sleep 1
        ls
        head $reference_genome
        determine_vaf.py $reference_genome $validation_bed $bam ${bam.simpleName}.snp.vcf ${bam.simpleName}.indel.vcf --threads ${task.cpus} 
        """
}

process FilterVariants{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple path(snp_vcf), path(indel_vcf), val(X), path(perbase_table)

    output:
        tuple path("${snp_vcf.simpleName}_filtered.snp.vcf"), path("${snp_vcf.simpleName}_filtered.indel.vcf")
    
    script:
        """
        vcf_filter.py -i $snp_vcf -o ${snp_vcf.simpleName}_filtered.snp.vcf -p $perbase_table \
        --min_dir_ratio $params.snp_filters.min_dir_ratio \
        --min_dir_count $params.snp_filters.min_dir_count \
        --min_dpq $params.snp_filters.min_dpq \
        --min_dpq_n $params.snp_filters.min_dpq_n \
        --min_dpq_ratio $params.snp_filters.min_dpq_ratio \
        --min_vaf $params.snp_filters.min_vaf \
        --min_rel_ratio $params.snp_filters.min_rel_ratio \
        --min_abq $params.snp_filters.min_abq

        vcf_filter.py -i $indel_vcf -o ${indel_vcf.simpleName}_filtered.indel.vcf -p $perbase_table \
        --min_dir_ratio $params.indel_filters.min_dir_ratio \
        --min_dir_count $params.indel_filters.min_dir_count \
        --min_dpq $params.indel_filters.min_dpq \
        --min_dpq_n $params.indel_filters.min_dpq_n \
        --min_dpq_ratio $params.indel_filters.min_dpq_ratio \
        --min_vaf $params.indel_filters.min_vaf \
        --min_rel_ratio $params.indel_filters.min_rel_ratio \
        --min_abq $params.indel_filters.min_abq
        """
}

process MergeNoisyVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple path(noisy_snp_vcf), path(noisy_indel_vcf)

    output:
        path("${noisy_snp_vcf.simpleName}.noisy_merged.vcf")
    
    script:
        """
	    mkdir tmpdir

        bcftools sort $noisy_snp_vcf -o ${noisy_snp_vcf.simpleName}.sorted.snp.vcf --temp-dir tmpdir
        bgzip ${noisy_snp_vcf.simpleName}.sorted.snp.vcf
        tabix ${noisy_snp_vcf.simpleName}.sorted.snp.vcf.gz

        bcftools sort $noisy_indel_vcf -o ${noisy_indel_vcf.simpleName}.sorted.indel.vcf --temp-dir tmpdir
        bgzip ${noisy_indel_vcf.simpleName}.sorted.indel.vcf
        tabix ${noisy_indel_vcf.simpleName}.sorted.indel.vcf.gz

        bcftools concat -a ${noisy_snp_vcf.simpleName}.sorted.snp.vcf.gz ${noisy_indel_vcf.simpleName}.sorted.indel.vcf.gz \
        -O v -o ${noisy_snp_vcf.simpleName}.noisy_merged.tmp.vcf
        bcftools sort ${noisy_snp_vcf.simpleName}.noisy_merged.tmp.vcf -o ${noisy_snp_vcf.simpleName}.noisy_merged.vcf --temp-dir tmpdir
        rm ${noisy_snp_vcf.simpleName}.noisy_merged.tmp.vcf

	    rm -r tmpdir
        """
}

process MergeFilteredVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple path(filtered_snp_vcf), path(filtered_indel_vcf)

    output:
        path("${filtered_snp_vcf.simpleName}.filtered_merged.vcf")
    
    script:
        """
	    mkdir tmpdir

        bcftools sort $filtered_snp_vcf -o ${filtered_snp_vcf.simpleName}.sorted.snp.vcf --temp-dir tmpdir
        bgzip ${filtered_snp_vcf.simpleName}.sorted.snp.vcf
        tabix ${filtered_snp_vcf.simpleName}.sorted.snp.vcf.gz

        bcftools sort $filtered_indel_vcf -o ${filtered_indel_vcf.simpleName}.sorted.indel.vcf --temp-dir tmpdir
        bgzip ${filtered_indel_vcf.simpleName}.sorted.indel.vcf
        tabix ${filtered_indel_vcf.simpleName}.sorted.indel.vcf.gz

        bcftools concat -a ${filtered_snp_vcf.simpleName}.sorted.snp.vcf.gz ${filtered_indel_vcf.simpleName}.sorted.indel.vcf.gz \
        -O v -o ${filtered_snp_vcf.simpleName}.filtered_merged.tmp.vcf
        bcftools sort ${filtered_snp_vcf.simpleName}.filtered_merged.tmp.vcf -o ${filtered_snp_vcf.simpleName}.filtered_merged.vcf --temp-dir tmpdir
        rm ${filtered_snp_vcf.simpleName}.filtered_merged.tmp.vcf

	    rm -r tmpdir
        """
}

process AnnotateVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        path(variant_vcf)


    output:
        path("${variant_vcf.simpleName}_annotated.vcf")
    
    script:
        """
        annotate_vcf.py $variant_vcf ${variant_vcf.simpleName}_annotated.vcf
        """
}
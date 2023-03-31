#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process minimap2 {
    queue 'middi'
    container "$params.azureRegistryServer/default/nwgs-minimap2:latest"
    cpus 32

    output:
        path 'sorted.bam'
        path 'sorted.bam.bai'

    script:
        """
        minimap2 --secondary=no --MD -ax map-ont -t 32 ${params.azureFileShare}/ref/${params.ref_genome} ${params.azureFileShare}/${params.reads} | samtools view -b -h -O "BAM" |  samtools sort -O "BAM" > sorted.bam
        samtools index sorted.bam
        """

    stub:
        """
        touch sorted.bam
        touch sorted.bam.bai
        """
}


process sniffles2 {
    queue 'middi'
    container "$params.azureRegistryServer/default/nwgs-sniffles2:latest"
    cpus 32
    
    input:
        path bam
        path index

    output:
        path 'sniffles.vcf'

    script:
        """
        sniffles --allow-overwrite --output-rnames -t 32 --minsvlen 10 --input $bam --vcf sniffles.vcf --reference ${params.azureFileShare}/ref/${params.ref_genome} --tandem-repeats ${params.azureFileShare}/ref/${params.tandem_repeat_bed}
        """

    stub:
        """
        touch sniffles.vcf
        """
}

process clair3 {
    queue 'middi'
    container "$params.azureRegistryServer/default/nwgs-clair3:latest"
    cpus 32

    input:
        path bam
        path index

    output:
        path 'phased_merge_output.vcf.gz'

    script:
        """
        run_clair3.sh --threads=32 \
        --include_all_ctgs \
        --bam_fn=$bam \
        --ref_fn=${params.azureFileShare}/ref/${params.ref_genome} \
        --platform=ont \
        --model_path=/root/miniconda3/envs/clair3/bin/models/r104_e81_sup_g5015 \
        --output=./ \
        --enable_phasing --longphase_for_phasing
        """

    stub:
        """
        touch phased_merge_output.vcf.gz
        """
}


process resultsout {
    queue 'middi'
    container "$params.azureRegistryServer/default/nwgs-bcftools:latest"
    cpus 8

    input:
        path sniffles2_vcf
        path clair3_vcf

    output:
        path 'sniffles.vcf.gz'
        path 'sniffles.vcf.gz.tbi'
        path 'minimap_on_target_clair_snvs.vcf.gz'
        path 'minimap_on_target_clair_snvs.vcf.gz.tbi'
        path 'minimap_on_target_clair_non-snvs.vcf.gz'
        path 'minimap_on_target_clair_non-snvs.vcf.gz.tbi'

    script:
        """
        cat <(grep ^# $sniffles2_vcf) <(grep -v ^# $sniffles2_vcf | LC_ALL=C sort -k1,1 -k2,2n) | sed '4i ##reference=hg38' | bgzip -c > sniffles.vcf.gz
        tabix -fp vcf sniffles.vcf.gz
        bcftools view --output-type z --types snps $clair3_vcf > minimap_on_target_clair_snvs.vcf.gz
        bcftools index -f --tbi minimap_on_target_clair_snvs.vcf.gz
        bcftools view --output-type z --exclude-types snps $clair3_vcf > nonsnv_tmp.vcf.gz
        bcftools norm --fasta-ref ${params.azureFileShare}/ref/${params.ref_genome} --output-type z ./nonsnv_tmp.vcf.gz > minimap_on_target_clair_non-snvs.vcf.gz
        bcftools index -f --tbi minimap_on_target_clair_non-snvs.vcf.gz
        """

    stub:
        """
        touch sniffles.vcf.gz
        """
}

process publishfiles {
    queue 'default'
    container "$params.azureRegistryServer/default/nwgs-bcftools:latest"
    cpus 2

    input:
        val sample_id
        path bam
        path bai
        path sniffles2_vcf
        path clair3_vcf
        path sniffles_gz
        path sniffles_gz_tbi
        path snv_vcf
        path snv_vcf_tbi
        path non_snv_vcf
        path non_snv_vcf_tbi

    script:
        """
        mkdir ${params.azureFileShare}/$sample_id
        cp $bam ${params.azureFileShare}/$sample_id/$sample_id-$bam.name
        cp $bai ${params.azureFileShare}/$sample_id/$sample_id-$bai.name
        cp $sniffles2_vcf ${params.azureFileShare}/$sample_id/$sample_id-$sniffles2_vcf.name
        cp $clair3_vcf ${params.azureFileShare}/$sample_id/$sample_id-$clair3_vcf.name
        cp $sniffles_gz ${params.azureFileShare}/$sample_id/$sample_id-$sniffles_gz.name
        cp $sniffles_gz_tbi ${params.azureFileShare}/$sample_id/$sample_id-$sniffles_gz_tbi.name
        cp $snv_vcf ${params.azureFileShare}/$sample_id/$sample_id-$snv_vcf.name
        cp $snv_vcf_tbi ${params.azureFileShare}/$sample_id/$sample_id-$snv_vcf_tbi.name
        cp $non_snv_vcf ${params.azureFileShare}/$sample_id/$sample_id-$non_snv_vcf.name
        cp $non_snv_vcf_tbi ${params.azureFileShare}/$sample_id/$sample_id-$non_snv_vcf_tbi.name
        """
}



workflow {
    sample_id = "$params.sample_id"

    minimap2()
    sniffles2(minimap2.out[0], minimap2.out[1])
    clair3(minimap2.out[0], minimap2.out[1])
    resultsout(sniffles2.out, clair3.out)
    publishfiles(sample_id, minimap2.out[0], minimap2.out[1], sniffles2.out, clair3.out, resultsout.out[0], resultsout.out[1], resultsout.out[2], resultsout.out[3], resultsout.out[4], resultsout.out[5])

}
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process minimap2 {
    queue 'default'
    container "$params.azureRegistryServer/default/nwgs-minimap2:latest"
    cpus 8

    publishDir "$params.azureFileShare/$params.outdir", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$filename" }

    input:
        val sample_id

    output:
        path 'sorted.bam'
        path 'sorted.bam.bai'

    script:
        """
        minimap2 --secondary=no --MD -ax map-ont -t 8 ${params.azureFileShare}/${params.ref_genome} ${params.azureFileShare}/${params.reads} | samtools view -b -h -O "BAM" |  samtools sort -O "BAM" > sorted.bam
        samtools index sorted.bam
        """

    stub:
        """
        touch sorted.bam
        touch sorted.bam.bai
        """
}


process sniffles2 {
    queue 'default'
    container "$params.azureRegistryServer/default/nwgs-sniffles2:latest"
    cpus 4

    publishDir "$params.azureFileShare/$params.outdir", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$filename" }
    
    input:
        val sample_id
        path bam
        path index

    output:
        path 'sniffles.vcf'

    script:
        """
        sniffles --allow-overwrite --output-rnames -t 4 --minsvlen 10 --input $bam --vcf sniffles.vcf --reference ${params.azureFileShare}/${params.ref_genome} --tandem-repeats ${params.azureFileShare}/${params.tandem_repeat_bed}
        """

    stub:
        """
        touch sniffles.vcf
        """
}

process clair3 {
    queue 'default'
    container "$params.azureRegistryServer/default/nwgs-clair3:latest"
    cpus 8

    publishDir "$params.azureFileShare/$params.outdir", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$filename" }

    input:
        val sample_id
        path bam
        path index

    output:
        path 'phased_merge_output.vcf.gz'

    script:
        """
        run_clair3.sh --threads=8 \
        --include_all_ctgs \
        --bam_fn=$bam \
        --ref_fn=${params.azureFileShare}/${params.ref_genome} \
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
    queue 'default'
    container "$params.azureRegistryServer/default/nwgs-bcftools:latest"
    cpus 4

    publishDir "$params.azureFileShare/$params.outdir", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$filename" }

    input:
        val sample_id
        path sniffles2_vcf
        path clair3_vcf

    output:
        path 'sniffles.vcf.gz*'
        path 'minimap_on_target_clair_snvs.vcf.gz*'
        path 'minimap_on_target_clair_non-snvs.vcf.gz*'

    script:
        """
        cat <(grep ^# $sniffles2_vcf) <(grep -v ^# $sniffles2_vcf | LC_ALL=C sort -k1,1 -k2,2n) | sed '4i ##reference=hg38' | bgzip -c > sniffles.vcf.gz
        tabix -fp vcf sniffles.vcf.gz
        bcftools view --output-type z --types snps $clair3_vcf > minimap_on_target_clair_snvs.vcf.gz
        bcftools index -f --tbi minimap_on_target_clair_snvs.vcf.gz
        bcftools view --output-type z --exclude-types snps $clair3_vcf > nonsnv_tmp.vcf.gz
        bcftools norm --fasta-ref ${params.azureFileShare}/${params.ref_genome} --output-type z ./nonsnv_tmp.vcf.gz > minimap_on_target_clair_non-snvs.vcf.gz
        bcftools index -f --tbi minimap_on_target_clair_non-snvs.vcf.gz
        """

    stub:
        """
        touch sniffles.vcf.gz
        """
}


workflow {
    sample_id = "$params.sample_id"

    minimap2(sample_id)
    sniffles2(sample_id, minimap2.out[0], minimap2.out[1])
    clair3(sample_id, minimap2.out[0], minimap2.out[1])
    resultsout(sample_id, sniffles2.out, clair3.out)

}
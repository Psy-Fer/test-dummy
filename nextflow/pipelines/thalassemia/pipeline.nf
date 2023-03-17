#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process minimap2 {
    container "$params.azureRegistryServer/default/nwgs-minimap2:latest"

    publishDir "$params.azureFileShare/$params.outdir", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$filename" }

    input:
        val sample_id
        val reads
        val ref

    output:
        path 'sorted.bam'
        path 'sorted.bam.bai'

    script:
        """
        minimap2 --secondary=no --MD -ax map-ont -t ${task.cpus} ${params.azureFileShare}/$ref ${params.azureFileShare}/$reads | samtools view -b -h -O "BAM" |  samtools sort -O "BAM" > sorted.bam
        samtools index sorted.bam
        """

    stub:
        """
        touch sorted.bam
        touch sorted.bam.bai
        """
}


process sniffles2 {
    container "$params.azureRegistryServer/default/nwgs-sniffles2:latest"

    publishDir "$params.azureFileShare/$params.outdir", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$filename" }
    
    input:
        val sample_id
        path bam
        path index
        val ref
        val trf

    output:
        path 'sniffles.vcf'

    script:
        """
        sniffles --allow-overwrite --output-rnames --minsvlen 10 --input $bam --vcf sniffles.vcf --reference ${params.azureFileShare}/$ref --tandem-repeats ${params.azureFileShare}/$trf
        """

    stub:
        """
        touch sniffles.vcf
        """
}

process clair3 {
    container "$params.azureRegistryServer/default/nwgs-clair3:latest"

    publishDir "$params.azureFileShare/$params.outdir", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$filename" }

    input:
        val sample_id
        path bam
        path index
        val ref
        val ref_index
        val trf

    output:
        path 'phased_merge_output.vcf.gz'

    script:
        """
        run_clair3.sh --threads=${task.cpus} \
        --include_all_ctgs \
        --bam_fn=$bam \
        --ref_fn=${params.azureFileShare}/$ref \
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
    container "$params.azureRegistryServer/default/nwgs-bcftools:latest"

    publishDir "$params.azureFileShare/$params.outdir", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id.$filename" }

    input:
        val sample_id
        path sniffles2_vcf
        path clair3_vcf
        val ref

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
        bcftools norm --fasta-ref ${params.azureFileShare}/$ref --output-type z ./nonsnv_tmp.vcf.gz > minimap_on_target_clair_non-snvs.vcf.gz
        bcftools index -f --tbi minimap_on_target_clair_non-snvs.vcf.gz
        """

    stub:
        """
        touch sniffles.vcf.gz
        """
}


workflow {
    reads = "$params.reads"
    ref = "$params.ref_genome"
    ref_index = "$params.ref_genome_index"
    trf = "$params.tandem_repeat_bed"
    sample_id = "$params.sample_id"

    minimap2(sample_id, reads, ref)
    sniffles2(sample_id, minimap2.out[0], minimap2.out[1], ref, trf)
    clair3(sample_id, minimap2.out[0], minimap2.out[1], ref, ref_index, trf)
    resultsout(sample_id, sniffles2.out, clair3.out, ref)

}
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process preprocess {
    queue 'default'
    container "$params.azureRegistryServer/default/nwgs-bcftools:latest"
    cpus 2

    output:
        path "$params.reads"
        path "$params.ref_genome"
        path "$params.ref_genome_index"
        path "$params.tandem_repeat_bed"
        path "$params.thal_regions"

    script:
        """
        echo "Staging files"
        ln -s ${params.azureFileShare}/${params.reads} ${params.reads}
        ln -s ${params.azureFileShare}/ref/${params.ref_genome} ${params.ref_genome}
        ln -s ${params.azureFileShare}/ref/${params.ref_genome_index} ${params.ref_genome_index}
        ln -s ${params.azureFileShare}/ref/${params.tandem_repeat_bed} ${params.tandem_repeat_bed}
        ln -s ${params.azureFileShare}/ref/${params.thal_regions} ${params.thal_regions}
        ls -la
        ls -la > ./ls_output.txt

        """

    stub:
        """
        touch ref.fna
        touch ref.fna.fai
        touch trf.bed
        touch regions.bed
        """
}


process minimap2 {
    queue 'middi'
    container "$params.azureRegistryServer/default/nwgs-minimap2:latest"
    cpus 32

    input:
        path reads
        path ref_genome
        path ref_genome_index

    output:
        path 'minimap2.sorted.bam'
        path 'minimap2.sorted.bam.bai'

    script:
        """
        minimap2 --secondary=no --MD -ax map-ont -t ${task.cpus} $ref_genome $reads | samtools sort - > minimap2.sorted.bam
        samtools index minimap2.sorted.bam
        """

    stub:
        """
        touch minimap2.sorted.bam
        touch minimap2.sorted.bam.bai
        """
}

process clair3 {
    queue 'middi'
    container "$params.azureRegistryServer/default/nwgs-clair3:latest"
    cpus 32

    input:
        val sample_id
        path ref_genome
        path ref_genome_index
        path bam
        path index
        path thal_regions

    output:
        path 'merge_output.vcf.gz'
        path 'merge_output.vcf.gz.tbi'

    script:
        """
        run_clair3.sh --threads=${task.cpus} \
        --sample_name=$sample_id \
        --ref_fn=$ref_genome \
        --bed_fn=$thal_regions \
        --bam_fn=$bam \
        --model_path=/root/miniconda3/envs/clair3/bin/models/r104_e81_sup_g5015 \
        --include_all_ctgs \
        --platform=ont \
        --output=./ \
        --snp_min_af=0.05 \
        --indel_min_af=0.10 \
        --min_mq=10 \
        --min_coverage=5

        """

    stub:
        """
        touch merge_output.vcf.gz
        touch merge_output.vcf.gz.tbi
        """
}

process whatshap {
    queue 'default'
    container "$params.azureRegistryServer/default/nwgs-clair3:latest"
    cpus 8

    input:
        val sample_id
        path ref_genome
        path ref_genome_index
        path bam
        path index
        path clair3_vcf
        path clair3_vcf_tbi

    output:
        path 'clair3_whatshap.vcf.gz'
        path 'clair3_whatshap.vcf.gz.tbi'
        path 'minimap_haplotagged.bam'
        path 'minimap_haplotagged.bam.bai'
        path 'clair3_whatshap.phase_stats.tsv'
        path 'clair3_whatshap.phase_blocks.gtf'
        path 'clair3_whatshap.haplotag_list.tsv'

    script:
        """
        whatshap phase \
        --sample $sample_id \
        --reference $ref_genome \
        --output clair3_whatshap.vcf.gz \
        --mapping-quality 10 \
        --indels \
        --ignore-read-groups \
        $clair3_vcf \
        $bam

        bcftools index --tbi clair3_whatshap.vcf.gz

        whatshap stats \
        --sample $sample_id \
        --tsv clair3_whatshap.phase_stats.tsv \
        --gtf clair3_whatshap.phase_blocks.gtf \
        clair3_whatshap.vcf.gz

        # Haplotag on-target alignments [whatshapp]
        
        whatshap haplotag \
        --sample $sample_id \
        --reference $ref_genome \
        --output minimap_haplotagged.bam \
        --tag-supplementary \
        --ignore-read-groups \
        --output-haplotag-list clair3_whatshap.haplotag_list.tsv \
        clair3_whatshap.vcf.gz \
        $bam
        
        samtools index minimap_haplotagged.bam

        """

    stub:
        """
        touch clair3_whatshap.vcf.gz
        touch clair3_whatshap.vcf.gz.tbi
        touch minimap_haplotagged.bam
        touch minimap_haplotagged.bam.bai
        touch clair3_whatshap.phase_stats.tsv
        touch clair3_whatshap.phase_blocks.gtf
        touch clair3_whatshap.haplotag_list.tsv
        """
}

process sniffles2 {
    queue 'middi'
    container "$params.azureRegistryServer/default/nwgs-sniffles2:latest"
    cpus 32
    
    input:
        val sample_id
        path ref_genome
        path ref_genome_index
        path tandem_repeat_bed
        path bam
        path index

    output:
        path 'sniffles.vcf.gz'
        path 'sniffles.snf'

    script:
        """
        sniffles --allow-overwrite \
        --threads ${task.cpus} \
        --sample-id $sample_id \
        --reference $ref_genome \
        --tandem-repeats $tandem_repeat_bed \
        --output-rnames \
        --minsvlen 20 \
        --input $bam \
        --vcf sniffles.vcf.gz \
        --snf sniffles.snf \
        --phase

        """

    stub:
        """
        touch sniffles.vcf.gz
        touch sniffles.snf
        """
}



process postproccess {
    queue 'middi'
    container "$params.azureRegistryServer/default/nwgs-clair3:latest"
    cpus 32

    input:
        path sniffles2_vcf
        path clair3_whatshap_vcf
        path clair3_whatshap_vcf_index
        path bam
        path index

    output:
        path 'sniffles.vcf.gz'
        path 'sniffles.vcf.gz.tbi'
        path 'clair3_whatshap_corrected.vcf.gz'
        path 'clair3_whatshap_corrected.vcf.gz.tbi'

    shell:
        """
        # index sniffles vcf
        bcftools index --tbi !{sniffles2_vcf}

        # correct zygosity of clair3 variants within sniffles deletions

        pypy3 /root/miniconda3/envs/clair3/bin/clair3.py SwitchZygosityBasedOnSVCalls \
        --threads !{task.cpus} \
        --bam_fn !{bam} \
        --clair3_vcf_input !{clair3_whatshap_vcf} \
        --sv_vcf_input !{sniffles2_vcf} \
        --vcf_output clair3_whatshap_corrected.vcf
        
        """

    stub:
        """
        touch sniffles.vcf.gz
        touch sniffles.vcf.gz.tbi
        touch clair3_whatshap_corrected.vcf.gz
        touch clair3_whatshap_corrected.vcf.gz.tbi
        """
}

process publishfiles {
    queue 'default'
    container "$params.azureRegistryServer/default/nwgs-bcftools:latest"
    cpus 2

    input:
        val sample_id
        path bam
        path bam_bai
        path clair3_whatshap_vcf
        path clair3_whatshap_vcf_tbi
        path minimap_haplotagged_bam
        path minimap_haplotagged_bam_bai
        path sniffles2_snf
        path clair3_whatshap_phase_stats
        path clair3_whatshap_phase_blocks_gtf
        path clair3_whatshap_haplotag_list
        path sniffles2_vcf
        path sniffles2_vcf_tbi
        path clair3_whatshap_corrected_vcf
        path clair3_whatshap_corrected_vcf_tbi

    script:
        """
        mkdir ${params.azureFileShare}/$sample_id
        cp $bam ${params.azureFileShare}/$sample_id/$sample_id-$bam.name
        cp $bam_bai ${params.azureFileShare}/$sample_id/$sample_id-$bam_bai.name
        cp $clair3_whatshap_vcf ${params.azureFileShare}/$sample_id/$sample_id-$clair3_whatshap_vcf.name
        cp $clair3_whatshap_vcf_tbi ${params.azureFileShare}/$sample_id/$sample_id-$clair3_whatshap_vcf_tbi.name
        cp $minimap_haplotagged_bam ${params.azureFileShare}/$sample_id/$sample_id-$minimap_haplotagged_bam.name
        cp $minimap_haplotagged_bam_bai ${params.azureFileShare}/$sample_id/$sample_id-$minimap_haplotagged_bam_bai.name
        cp $sniffles2_snf ${params.azureFileShare}/$sample_id/$sample_id-$sniffles2_snf.name
        cp $clair3_whatshap_phase_stats ${params.azureFileShare}/$sample_id/$sample_id-$clair3_whatshap_phase_stats.name
        cp $clair3_whatshap_phase_blocks_gtf ${params.azureFileShare}/$sample_id/$sample_id-$clair3_whatshap_phase_blocks_gtf.name
        cp $clair3_whatshap_haplotag_list ${params.azureFileShare}/$sample_id/$sample_id-$clair3_whatshap_haplotag_list.name
        cp $sniffles2_vcf ${params.azureFileShare}/$sample_id/$sample_id-$sniffles2_vcf.name
        cp $sniffles2_vcf_tbi ${params.azureFileShare}/$sample_id/$sample_id-$sniffles2_vcf_tbi.name
        cp $clair3_whatshap_corrected_vcf ${params.azureFileShare}/$sample_id/$sample_id-$clair3_whatshap_corrected_vcf.name
        cp $clair3_whatshap_corrected_vcf_tbi ${params.azureFileShare}/$sample_id/$sample_id-$clair3_whatshap_corrected_vcf_tbi.name
        """

    stub:
        """
        touch minimap2.sorted.bam
        touch minimap2.sorted.bam.bai
        touch clair3_whatshap.vcf.gz
        touch clair3_whatshap.vcf.gz.tbi
        touch minimap_haplotagged.bam
        touch minimap_haplotagged.bam.bai
        touch sniffles.snf
        touch clair3_whatshap.phase_stats.tsv
        touch clair3_whatshap.phase_blocks.gtf
        touch clair3_whatshap.haplotag_list.tsv
        touch sniffles.vcf.gz
        touch sniffles.vcf.gz.tbi
        touch clair3_whatshap_corrected.vcf.gz
        touch clair3_whatshap_corrected.vcf.gz.tbi
        """
}



workflow {

    sample_id = "$params.sample_id"

    preprocess()
    minimap2(preprocess.out[0], preprocess.out[1], preprocess.out[2])
    clair3(sample_id,  preprocess.out[1], preprocess.out[2], minimap2.out[0], minimap2.out[1], preprocess.out[4])
    whatshap(sample_id, preprocess.out[1], preprocess.out[2], minimap2.out[0], minimap2.out[1], clair3.out[0], clair3.out[1])
    sniffles2(sample_id, preprocess.out[1], preprocess.out[2], preprocess.out[3], whatshap.out[2], whatshap.out[3])
    postproccess(sniffles2.out[0], whatshap.out[0], whatshap.out[1], whatshap.out[2], whatshap.out[3])
    publishfiles(sample_id, minimap2.out[0], minimap2.out[1], sniffles2.out[1], whatshap.out[0], whatshap.out[1], whatshap.out[2],
                 whatshap.out[3], whatshap.out[4], whatshap.out[5], whatshap.out[6], postproccess.out[0], 
                 postproccess.out[1], postproccess.out[2], postproccess.out[3])

}

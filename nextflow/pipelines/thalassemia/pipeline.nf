#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process preprocess {
    queue 'default'
    container "$params.azureRegistryServer/default/nwgs-bcftools:latest"
    cpus 2

    input:
        path reads
        path ref
        path ref_index
        path trf
        path bed


    output:
        path "$reads.name"
        path "$ref.name"
        path "$ref_index.name"
        path "$trf.name"
        path "$bed.name"

    script:
        """
        echo "Staging files"

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
        path thal_region
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
        path thal_region
        path bam
        path index

    output:
        path 'sniffles.vcf.gz'
        path 'sniffles.snf'

    script:
        """
        sniffles --allow-overwrite \
        -threads ${task.cpus} \
        --sample-id $sample_id \
        --reference $ref_genome \
        --tandem-repeats $tandem_repeat_bed \
        --output-rnames \
        --minsvlen 20 \
        --input $bam \
        --vcf sniffles.vcf \
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
        path 'combined_variants.vcf.gz'
        path 'combined_variants.vcf.gz.tbi'

    shell:
        """

        # bcftools view --header-only !{clair3_whatshap_vcf} > tmp.vcf
        # insertions
        # bcftools view --no-header --output-type v --types indels !{clair3_whatshap_vcf} | awk 'length(\$4) == 1 && length(\$5) > 1 && length(\$5) < 21' >> tmp.vcf
        # deletions
        # bcftools view --no-header --output-type v --types indels !{clair3_whatshap_vcf} | awk 'length(\$5) == 1 && length(\$4) > 1 && length(\$4) < 21' >> tmp.vcf
        # variants that are not indels
        # bcftools view --no-header --output-type v --exclude-types indels !{clair3_whatshap_vcf} >> tmp.vcf
        # sort
        # bcftools sort --output-type z tmp.vcf > clair3_whatshap_filtered.vcf.gz
        # index
        # bcftools index --tbi clair3_whatshap_filtered.vcf.gz

        # index sniffles vcf
        bcftools index --tbi !{sniffles2_vcf}

        # correct zygosity of clair3 variants within sniffles deletions

        pypy3 /root/miniconda3/envs/clair3/bin/clair3.py SwitchZygosityBasedOnSVCalls \
        --threads !{task.cpus} \
        --bam_fn !{bam} \
        --clair3_vcf_input !{clair3_whatshap_vcf} \
        --sv_vcf_input !{sniffles2_vcf} \
        --vcf_output clair3_whatshap_corrected.vcf

        # bcftools index --tbi clair3_whatshap_corrected.vcf.gz
            
        # Combine SVs and small variants into a single phased VCF [bcftools merge]
        # different formatting of sniffles and clair3 VCFs might cause some problems here: maybe not possible to just merge them straight up.

        bcftools merge \
        --force-samples \
        --merge both \
        --output-type z \
        --output ./combined_variants.vcf.gz \
        clair3_whatshap_corrected.vcf.gz \
        !{sniffles2_vcf}

        bcftools index --tbi combined_variants.vcf.gz

        """

    stub:
        """
        touch sniffles.vcf.gz
        touch sniffles.vcf.gz.tbi
        # touch clair3_whatshap_filtered.vcf.gz
        # touch clair3_whatshap_filtered.vcf.gz.tbi
        touch clair3_whatshap_corrected.vcf.gz
        touch clair3_whatshap_corrected.vcf.gz.tbi
        touch combined_variants.vcf.gz
        touch combined_variants.vcf.gz.tbi
        """
}

process publishfiles {
    queue 'default'
    container "$params.azureRegistryServer/default/nwgs-bcftools:latest"
    cpus 2

    publishDir "$params.azureFileShare/$params.sample_id/pubdir", mode: 'copy', overwrite: true, saveAs: { filename -> "$sample_id-$filename" }

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
        path combined_variants_vcf
        path combined_variants_vcf_tbi
    
    output:
        path 'minimap2.sorted.bam'
        path 'minimap2.sorted.bam.bai'
        path 'clair3_whatshap.vcf.gz'
        path 'clair3_whatshap.vcf.gz.tbi'
        path 'minimap_haplotagged.bam'
        path 'minimap_haplotagged.bam.bai'
        path 'sniffles.snf'
        path 'clair3_whatshap.phase_stats.tsv'
        path 'clair3_whatshap.phase_blocks.gtf'
        path 'clair3_whatshap.haplotag_list.tsv'
        path 'sniffles.vcf.gz'
        path 'sniffles.vcf.gz.tbi'
        path 'clair3_whatshap_filtered.vcf.gz'
        path 'clair3_whatshap_filtered.vcf.gz.tbi'
        path 'clair3_whatshap_corrected.vcf.gz'
        path 'clair3_whatshap_corrected.vcf.gz.tbi'
        path 'combined_variants.vcf.gz'
        path 'combined_variants.vcf.gz.tbi'

    script:
        """
        mkdir ${params.azureFileShare}}/$sample_id
        cp $bam ${params.azureFileShare}}/$sample_id/$sample_id-$bam.name
        cp $bam_bai ${params.azureFileShare}}/$sample_id/$sample_id-$bam_bai.name
        cp $clair3_whatshap_vcf ${params.azureFileShare}}/$sample_id/$sample_id-$clair3_whatshap_vcf.name
        cp $clair3_whatshap_vcf_tbi ${params.azureFileShare}}/$sample_id/$sample_id-$clair3_whatshap_vcf_tbi.name
        cp $minimap_haplotagged_bam ${params.outdir}/$sample_id/$sample_id-$minimap_haplotagged_bam.name
        cp $minimap_haplotagged_bam_bai ${params.azureFileShare}}/$sample_id/$sample_id-$minimap_haplotagged_bam_bai.name
        cp $sniffles2_snf ${params.azureFileShare}}/$sample_id/$sample_id-$sniffles2_snf.name
        cp $clair3_whatshap_phase_stats ${params.azureFileShare}}/$sample_id/$sample_id-$clair3_whatshap_phase_stats.name
        cp $clair3_whatshap_phase_blocks_gtf ${params.azureFileShare}}/$sample_id/$sample_id-$clair3_whatshap_phase_blocks_gtf.name
        cp $clair3_whatshap_haplotag_list ${params.azureFileShare}}/$sample_id/$sample_id-$clair3_whatshap_haplotag_list.name
        cp $sniffles2_vcf ${params.azureFileShare}}/$sample_id/$sample_id-$sniffles2_vcf.name
        cp $sniffles2_vcf_tbi ${params.azureFileShare}}/$sample_id/$sample_id-$sniffles2_vcf_tbi.name
        cp $clair3_whatshap_corrected_vcf ${params.azureFileShare}}/$sample_id/$sample_id-$clair3_whatshap_corrected_vcf.name
        cp $clair3_whatshap_corrected_vcf_tbi ${params.azureFileShare}}/$sample_id/$sample_id-$clair3_whatshap_corrected_vcf_tbi.name
        cp $combined_variants_vcf ${params.azureFileShare}}/$sample_id/$sample_id-$combined_variants_vcf.name
        cp $combined_variants_vcf_tbi ${params.azureFileShare}}/$sample_id/$sample_id-$combined_variants_vcf_tbi.name
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
        touch combined_variants.vcf.gz
        touch combined_variants.vcf.gz.tbi
        """
}



workflow {

    sample_id = "$params.sample_id"

    reads = file("$params.azureFileShare/ref/$params.reads")
    ref_genome = file("$params.azureFileShare/ref/$params.ref_genome")
    ref_genome_index = file("$params.azureFileShare/ref/$params.ref_genome_index")
    tandem_repeat_bed = file("$params.azureFileShare/ref/$params.tandem_repeat_bed")
    thal_regions = file("$params.azureFileShare/ref/$params.thal_regions")

    preprocess(reads, ref_genome, ref_genome_index, tandem_repeat_bed, thal_regions)
    minimap2(preprocess.out[0], preprocess.out[1], preprocess.out[2])
    clair3(sample_id,  preprocess.out[1], preprocess.out[2], minimap2.out[0], minimap2.out[1], preprocess.out[4])
    whatshap(sample_id, preprocess.out[1], preprocess.out[2], minimap2.out[0], minimap2.out[1], preprocess.out[4], clair3.out[0], clair3.out[1])
    sniffles2(sample_id, preprocess.out[1], preprocess.out[2], preprocess.out[3], whatshap.out[2], whatshap.out[3])
    postproccess(sniffles2.out[0], whatshap.out[0], whatshap.out[1], whatshap.out[2], whatshap.out[3])
    publishfiles(sample_id, minimap2.out[0], minimap2.out[1], sniffles2.out[1], whatshap.out[0], whatshap.out[1], whatshap.out[2],
                 whatshap.out[3], whatshap.out[4], whatshap.out[5], whatshap.out[6], postproccess.out[0], 
                 postproccess.out[1], postproccess.out[2], postproccess.out[3], postproccess.out[4], postproccess.out[5],
                 postproccess.out[6], postproccess.out[7])

}

#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Utils modules
include { AddReplaceReadgroup as Samtools_AddReplaceReadgroup } from './Modules/Samtools/1.15/AddReplaceReadgroup.nf'
include { Annotate as Bedtools_Annotate_Clair3 } from './Modules/bedtools/1.15.1--h0ea216a_0/Annotate.nf'
include { Annotate as Bedtools_Annotate_Bed } from './Modules/bedtools/1.15.1--h0ea216a_0/Annotate.nf'
include { Annotate as Bedtools_Annotate_Region } from './Modules/bedtools/1.15.1--h0ea216a_0/Annotate.nf'
include { CollectMultipleMetrics as PICARD_CollectMultipleMetrics } from './Modules/Picard/2.26.4/CollectMultipleMetrics.nf' params(genome:"$params.genome_fasta", optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE")
include { CollectWgsMetrics as PICARD_CollectWgsMetrics } from './Modules/Picard/2.26.4/CollectWgsMetrics.nf' params(genome:"$params.genome_fasta", optional: "MINIMUM_MAPPING_QUALITY=1 MINIMUM_BASE_QUALITY=1 ")
include { ExportParams as Workflow_ExportParams } from './Modules/Utils/workflow.nf'
include { Fastq as Samtools_Fastq } from './Modules/Samtools/1.15/Fastq.nf' params(tags: " -T RG,Mm,Ml,MM,ML ", , roi: params.roi)
include { FilterCondition as Sambamba_Filter_Condition } from './Modules/Sambamba/1.0.0/Filter.nf' params(conditions: params.conditions)
include { FilterHaplotypePhaseset as Sambamba_Filter_Haplotype_Phaseset } from './Modules/Sambamba/1.0.0/Filter.nf'
include { FilterPairs as Duplex_FilterPairs } from './Modules/duplex_tools/0.2.17/FilterPairs.nf'
include { FilterSamReads as PICARD_FilterSamReads } from './Modules/Picard/2.26.4/FilterSamReads.nf' params(optional: " FILTER=excludeReadList")
include { FilterVcfs as GATK_FilterSNV_Target_Bed } from './Modules/GATK/4.2.1.0/FilterVCFs.nf' params(genome: params.genome_fasta, filter: "SNP")
include { FilterVcfs as GATK_FilterSNV_Target_Region } from './Modules/GATK/4.2.1.0/FilterVCFs.nf' params(genome: params.genome_fasta, filter: "SNP")
include { GetPhaseSet } from './Modules/GetPhaseSet.nf'
include { HaplotypeCaller_SMN as GATK_HaplotypeCaller_Bed } from './Modules/GATK/4.2.1.0/HaplotypeCaller.nf' params(genome: params.genome_fasta, compress: true, extention: "_bed", optional:"--intervals $params.calling_target_bed --dont-use-soft-clipped-bases --pair-hmm-implementation  LOGLESS_CACHING")
include { HaplotypeCaller_SMN as GATK_HaplotypeCaller_Region } from './Modules/GATK/4.2.1.0/HaplotypeCaller.nf' params(genome: params.genome_fasta, compress: true, extention: "_region", optional:"--intervals $params.calling_target_region --dont-use-soft-clipped-bases --pair-hmm-implementation  LOGLESS_CACHING")
include { Haplotag as Whatshap_Haplotag_Target_Bed } from './Modules/Whatshap/1.7/Haplotag.nf' params (genome: params.genome_fasta, extention: "_bed")
include { Haplotag as Whatshap_Haplotag_Target_Region } from './Modules/Whatshap/1.7/Haplotag.nf' params (genome: params.genome_fasta, extention: "_region")
include { Index as Sambamba_Index_Longshot } from './Modules/Sambamba/1.0.0/Index.nf'
include { Index as Sambamba_Index_Deduplex } from './Modules/Sambamba/1.0.0/Index.nf'
include { Index as Sambamba_Index_ReadGroup } from './Modules/Sambamba/1.0.0/Index.nf'
include { Index as Sambamba_Index_Target_Bed } from './Modules/Sambamba/1.0.0/Index.nf'
include { Index as Sambamba_Index_Target_Region } from './Modules/Sambamba/1.0.0/Index.nf'
include { Index as Sambamba_Index_Merge } from './Modules/Sambamba/1.0.0/Index.nf'
include { LongshotPhase } from './Modules/Longshot/0.4.1/Phase.nf'
include { Mapping as Minimap2_remap } from './Modules/Minimap2/2.26--he4a0461_1/Mapping.nf' params(optional: " -y -ax map-ont", genome_fasta: params.genome_fasta)
include { Merge as Samtools_Merge } from './Modules/Samtools/1.15/Merge.nf'
include { MultiQC } from './Modules/MultiQC/1.10/MultiQC.nf' params(optional: "--config $baseDir/assets/multiqc_config.yaml")
include { PairsFromSummary as Duplex_PairsFromSummary } from './Modules/duplex_tools/0.2.17/PairsFromSummary.nf'
include { Phase as Whatshap_Phase_Target_Bed } from './Modules/Whatshap/1.7/Phase.nf' params (genome: params.genome_fasta)
include { Phase as Whatshap_Phase_Target_Region } from './Modules/Whatshap/1.7/Phase.nf' params (genome: params.genome_fasta)
include { ReBasecallingGuppy } from './Modules/GuppyBasecalling.nf'
include { ViewSort as Sambamba_ViewSort_remap } from './Modules/Sambamba/1.0.0/ViewSort.nf'
include { VariantCaller as Clair3_VariantCaller } from './Modules/Clair3/1.0.4--py39hf5e1c6e_3/VariantCaller.nf' params(
    genome: "$params.genome_fasta",
    clair3model: "$params.clair3model",
    optional: " --haploid_precise --platform=ont --enable_long_indel"
)
include { VariantCaller as Sniffles2_VariantCaller } from './Modules/Sniffles2/2.2--pyhdfd78af_0/VariantCaller.nf' params(optional: "")
include { VariantFiltrationSnpIndel as GATK_VariantFiltration_Clair3 } from './Modules/GATK/4.2.1.0/VariantFiltration.nf' params(
    genome: "$params.genome_fasta", snp_filter: "$params.clair3_snp_filter",
    snp_cluster: "$params.clair3_snp_cluster", indel_filter: "$params.gatk_indel_filter", compress: true
)
include { VariantFiltrationSnpIndel as GATK_VariantFiltration_Bed } from './Modules/GATK/4.2.1.0/VariantFiltration.nf' params(
    genome: "$params.genome_fasta", snp_filter: "$params.gatk_snp_filter",
    snp_cluster: "$params.gatk_snp_cluster", indel_filter: "$params.gatk_indel_filter", compress: true
)
include { VariantFiltrationSnpIndel as GATK_VariantFiltration_Region } from './Modules/GATK/4.2.1.0/VariantFiltration.nf' params(
    genome: "$params.genome_fasta", snp_filter: "$params.gatk_snp_filter",
    snp_cluster: "$params.gatk_snp_cluster", indel_filter: "$params.gatk_indel_filter", compress: true
)
include { VersionLog } from './Modules/VersionLog.nf'
include { ZipIndex as Tabix_Zip_Index_Bed } from './Modules/Tabix/1.11/Index.nf'
include { ZipIndex as Tabix_Zip_Index_Region } from './Modules/Tabix/1.11/Index.nf'
include { ZipIndex as Tabix_Zip_Index_Bedtools_Clair3 } from './Modules/Tabix/1.11/Index.nf'
include { ZipIndex as Tabix_Zip_Index_Bedtools_Bed } from './Modules/Tabix/1.11/Index.nf'
include { ZipIndex as Tabix_Zip_Index_Bedtools_Region } from './Modules/Tabix/1.11/Index.nf'


def analysis_id = params.outdir.split('/')[-1]
sample_id = params.sample_id


if (params.method == "SMA_adaptive"){
    ploidy_list = Channel.of(1..params.ploidy)
}

workflow {
    if( params.start == 'bam' ){
        // Get fast5 and mapped bams from input folder
        fast5_files = Channel.fromPath(params.input_path +  "/workspace/fast5_pass/*.fast5").toList()
        bam_files = Channel.fromPath(params.input_path +  "/pass/*.bam").toList()
        summary_file = Channel.fromPath(params.input_path +  "/sequencing_summary.txt").toList()
    }
    else if( params.start == 'bam_single' || params.start == 'bam_single_remap'){
        //Get BAM file, and only BAM file as fast5 and summary are not available
        bam_file = Channel.fromPath(params.input_path +  "/*.bam").toList()
    }
    else if( params.start == 'rebase' ){
        //Re-basecalling
        fast5 = Channel.fromPath(params.input_path +  "/fast5_*/*.fast5").toList()
        ReBasecallingGuppy(params.input_path, sample_id, fast5.sum{it.size()})
        fast5_files = ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files, summary_file -> fast5_files}
        bam_files = ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files, summary_file -> bam_files}
        summary_file = ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files, summary_file -> summary_file}
    }
    else if( params.start == 'bam_remap' ){
        // Get fast5 and mapped bams from input folder
        fast5_files = Channel.fromPath(params.input_path +  "/workspace/fast5_pass/*.fast5").toList()
        bam_files = Channel.fromPath(params.input_path +  "/pass/*.bam").toList()
        summary_file = Channel.fromPath(params.input_path +  "/sequencing_summary.txt").toList()
    }
    else{
        error  """
            Invalid alignment mode: ${params.start}.
            This should be bam (start from basecalled data),
            bam_remap (start from bam, but perform remapping with minimap2),
            bam_single (start from single bam without sequencing information),
            bam_single_remap (start from single bam without sequencing information and perform remapping),
            or rebase (full re-basecalling)
        """
    }


    if( params.start == 'bam_single' ||  params.start == 'bam_single_remap' ){
        // Index MergeSort BAM
        Sambamba_Index_Merge(bam_file.map{bam_file -> [sample_id, bam_file]})

        // Filter for minimum readlength
        Sambamba_Filter_Condition(bam_file.combine(Sambamba_Index_Merge.out.map{sample_id, bai_file -> bai_file}))

        bam_file = Sambamba_Filter_Condition.out
    }
    else{
        // MergeSort BAMs
        Samtools_Merge(bam_files.map{bam_files -> [sample_id, bam_files]})

        // Index MergeSort BAM
        Sambamba_Index_Merge(Samtools_Merge.out)

        // Filter for minimum readlength
        Sambamba_Filter_Condition(Samtools_Merge.out
            .map{ sample_id, bam_file -> bam_file }
            .combine(Sambamba_Index_Merge.out
                .map{ sample_id, bai_file -> bai_file }
            )
        )

        // Identify readpairs
        Duplex_PairsFromSummary(sample_id, summary_file)

        // Identify possible duplex reads from read pairs
        Duplex_FilterPairs(Duplex_PairsFromSummary.out, Sambamba_Filter_Condition.out)

        //Filter BAM for duplicate duplex read
        PICARD_FilterSamReads(Sambamba_Filter_Condition.out, Duplex_FilterPairs.out)

        //Index BAM file
        Sambamba_Index_Deduplex(PICARD_FilterSamReads.out.map{bam_file ->[sample_id, bam_file]})

        bam_file = PICARD_FilterSamReads.out.combine(
            Sambamba_Index_Deduplex.out.map{sample_id, bai_file -> bai_file}
        )

    }


    if( params.start == 'bam_remap' || params.start == 'bam_single_remap' ){
        // Extract FASTQ from BAM
        Samtools_Fastq(bam_file)

         // Re-map ROI fastq
        Minimap2_remap(Samtools_Fastq.out)

        // Sort SAM to BAM
        Sambamba_ViewSort_remap(Minimap2_remap.out.map{fastq, sam_file -> [sample_id, fastq , sam_file]})

        bam_file = Sambamba_ViewSort_remap.out.map{sample_id, rg_id, bam_file, bai_file -> [bam_file, bai_file]}

    }


    // Add readgroup to BAM
    Samtools_AddReplaceReadgroup(sample_id, bam_file)

    // Index readgroup BAM
    Sambamba_Index_ReadGroup(Samtools_AddReplaceReadgroup.out.map{bam_file -> [sample_id, bam_file]})

    bam_file = Samtools_AddReplaceReadgroup.out.combine(
        Sambamba_Index_ReadGroup.out.map{sample_id, bai_file -> bai_file})
        .map{bam_file, bai_file -> [sample_id, bam_file, bai_file]}


    if (params.method == "wgs"){
        //Phasing BAM
        LongshotPhase(bam_file.map{sample_id, bam_file, bai_file -> [bam_file, bai_file]})

        // BAMIndex
        Sambamba_Index_Longshot(LongshotPhase.out.map{bam_file, vcf_file -> [sample_id, bam_file]})

    }


    if (params.method == "SMA_adaptive"){
        // Variant calling
        GATK_HaplotypeCaller_Bed(
            bam_file.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file, params.ploidy]}
        )

        // Filter SNV only Bed-file variants
        GATK_FilterSNV_Target_Bed(GATK_HaplotypeCaller_Bed.out)

        // Whatshapp polyphase Bed-file variants
        Whatshap_Phase_Target_Bed(GATK_FilterSNV_Target_Bed.out)

        // bgzip and index VCF Bed-file variants
        Tabix_Zip_Index_Bed(
            Whatshap_Phase_Target_Bed.out.map{sample_id, vcf_file, ploidy -> [sample_id, vcf_file]}
        )

        //Annotate Homopolymer VCF
        Bedtools_Annotate_Bed(
            Tabix_Zip_Index_Bed.out.map{sample_id, vcf_file, vcf_file_index -> [vcf_file, vcf_file_index]}
        )

        //Index Homopolymer annotated VCF file
        Tabix_Zip_Index_Bedtools_Bed(Bedtools_Annotate_Bed.out.map{vcf_file -> [sample_id, vcf_file]})

        //Filter VCF
        GATK_VariantFiltration_Bed(Tabix_Zip_Index_Bedtools_Bed.out)

        // Whatshapp haplotag Bed-file variants
        Whatshap_Haplotag_Target_Bed(
            GATK_HaplotypeCaller_Bed.out
            .map{sample_id, bam_file, bai_file, vcf_file, vcf_index, ploidy -> [sample_id, bam_file, bai_file, ploidy]}
            .join(GATK_VariantFiltration_Bed.out)
        )

        // Index BAM file and publish Bed-file variants
        Sambamba_Index_Target_Bed(Whatshap_Haplotag_Target_Bed.out.map{bam_file -> [sample_id, bam_file]})

        // Variant calling on used defined region
        GATK_HaplotypeCaller_Region(
            bam_file.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file, params.ploidy]}
        )

        //Filter SNV only region variants
        GATK_FilterSNV_Target_Region(GATK_HaplotypeCaller_Region.out)

        // Whatshapp polyphase region
        Whatshap_Phase_Target_Region(GATK_FilterSNV_Target_Region.out)

        // bgzip and index VCF
        Tabix_Zip_Index_Region(
            Whatshap_Phase_Target_Region.out.map{sample_id, vcf_file, ploidy -> [sample_id, vcf_file]}
        )

        //Annotate Homopolymer VCF
        Bedtools_Annotate_Region(
            Tabix_Zip_Index_Region.out.map{sample_id, vcf_file, vcf_file_index -> [vcf_file, vcf_file_index]}
        )

        //Index Homopolymer annotated VCF file
        Tabix_Zip_Index_Bedtools_Region(Bedtools_Annotate_Region.out.map{vcf_file -> [sample_id, vcf_file]})

        //Filter VCF
        GATK_VariantFiltration_Region(Tabix_Zip_Index_Bedtools_Region.out)

        // Whatshapp haplotag
        Whatshap_Haplotag_Target_Region(
            GATK_HaplotypeCaller_Region.out
            .map{sample_id, bam_file, bai_file, vcf_file, vcf_index, ploidy -> [sample_id, bam_file, bai_file, ploidy]}
            .join(GATK_VariantFiltration_Region.out)
        )

        // Index BAM file and publish region variants
        Sambamba_Index_Target_Region(Whatshap_Haplotag_Target_Region.out.map{bam_file -> [sample_id, bam_file]})


        // Get correct phasegroup
        GetPhaseSet(GATK_VariantFiltration_Bed.out)

        // Split BAM into single haplotypes
        Sambamba_Filter_Haplotype_Phaseset(Whatshap_Haplotag_Target_Bed.out.combine(
            Sambamba_Index_Target_Bed.out.map{sample_id, bai_file -> bai_file}
        ).combine(
            ploidy_list
        ).combine(
           GetPhaseSet.out
        ))

        //Clair3 calling on haplotype BAMs
        Clair3_VariantCaller(Sambamba_Filter_Haplotype_Phaseset.out)

        // Annotate Clair3 VCF
        Bedtools_Annotate_Clair3(Clair3_VariantCaller.out)

        //Index Homopolymer annotated VCF file
        Tabix_Zip_Index_Bedtools_Clair3(Bedtools_Annotate_Clair3.out.map{vcf_file -> [sample_id, vcf_file]})

        //Filter VCF
        GATK_VariantFiltration_Clair3(Tabix_Zip_Index_Bedtools_Clair3.out)

        //Sniffles SV variant calling on haplotype BAMs
        Sniffles2_VariantCaller(Sambamba_Filter_Haplotype_Phaseset.out)

    }

    // QC stats
    PICARD_CollectMultipleMetrics(bam_file)
    PICARD_CollectWgsMetrics(bam_file)
    MultiQC(analysis_id, Channel.empty().mix(
        PICARD_CollectMultipleMetrics.out,
        PICARD_CollectWgsMetrics.out
    ).collect())

    // Create log files: Repository versions and Workflow params
    VersionLog()
    Workflow_ExportParams()
}

// Workflow completion notification
workflow.onComplete {
    // HTML Template
    def template = new File("$baseDir/assets/workflow_complete.html")
    def binding = [
        runName: analysis_id,
        workflow: workflow
    ]
    def engine = new groovy.text.GStringTemplateEngine()
    def email_html = engine.createTemplate(template).make(binding).toString()

    // Send email
    if (workflow.success) {
        def subject = "ONT Workflow Successful: ${analysis_id}"
        sendMail(
            to: params.email.trim(),
            subject: subject,
            body: email_html,
            attach: "${params.outdir}/QC/${analysis_id}_multiqc_report.html"
        )

    } else {
        def subject = "ONT Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}


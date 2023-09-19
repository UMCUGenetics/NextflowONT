#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Utils modules
include ExportParams as Workflow_ExportParams from './NextflowModules/Utils/workflow.nf'

//include Index as Sambamba_Index from './NextflowModules/Sambamba/0.7.0/Index.nf'
include AddReadgroup as Samtools_AddReadgroup from './NextflowModules/Samtools/1.15/AddReadgroup.nf'
include Merge as Sambamba_Merge from './NextflowModules/Sambamba/0.7.0/Merge.nf' 
include Filter_ROI as Sambamba_Filter_ROI from './NextflowModules/Sambamba/0.7.0/Filter.nf'
include Filter_Condition as Sambamba_Filter_Condition from './NextflowModules/Sambamba/0.7.0/Filter.nf' params(conditions: params.conditions)
include Split as Sambamba_Split from './NextflowModules/Sambamba/0.7.0/Split.nf'
include Index as Sambamba_Index_Longshot from './NextflowModules/Sambamba/0.7.0/Index.nf'
include Index as Sambamba_Index_Deduplex from './NextflowModules/Sambamba/0.7.0/Index.nf'
include LongshotPhase from './NextflowModules/Longshot/0.4.1/Phase.nf'
include PairsFromSummary as Duplex_PairsFromSummary from './NextflowModules/duplex_tools/0.2.17/PairsFromSummary.nf'
include FilterPairs as Duplex_FilterPairs from './NextflowModules/duplex_tools/0.2.17/FilterPairs.nf'
include FilterSamReads as PICARD_FilterSamReads from './NextflowModules/Picard/2.26.4/FilterSamReads.nf' params(optional: " FILTER=excludeReadList")

include ToSAM as Sambamba_ToSam from './NextflowModules/Sambamba/0.7.0/ToSAM.nf'
include ToSAM as Sambamba_ToSam_hap1 from './NextflowModules/Sambamba/0.7.0/ToSAM.nf'
include ToSAM as Sambamba_ToSam_hap2 from './NextflowModules/Sambamba/0.7.0/ToSAM.nf'
include ToSAM as Sambamba_ToSam_nohap from './NextflowModules/Sambamba/0.7.0/ToSAM.nf'
include GetReadIDs as Sambamba_GetReadIDs from './NextflowModules/Sambamba/0.7.0/GetReads.nf'
include GetReadIDs as Sambamba_GetReadIDs_hap1 from './NextflowModules/Sambamba/0.7.0/GetReads.nf'
include GetReadIDs as Sambamba_GetReadIDs_hap2 from './NextflowModules/Sambamba/0.7.0/GetReads.nf'
include GetReadIDs as Sambamba_GetReadIDs_nohap from './NextflowModules/Sambamba/0.7.0/GetReads.nf'
include Index as STRiqueIndex from './NextflowModules/STRique/0.4.2/Index.nf'
include CallRepeat as STRiqueCallRepeat from './NextflowModules/STRique/0.4.2/CallRepeats.nf'
include CallRepeat as STRiqueCallRepeat_hap1 from './NextflowModules/STRique/0.4.2/CallRepeats.nf'
include CallRepeat as STRiqueCallRepeat_hap2 from './NextflowModules/STRique/0.4.2/CallRepeats.nf'
include CallRepeat as STRiqueCallRepeat_nohap from './NextflowModules/STRique/0.4.2/CallRepeats.nf'

include Fastq as Samtools_Fastq from './NextflowModules/Samtools/1.15/Fastq.nf' params(tags: " -T RG,Mm,Ml ", , roi: params.roi)
include Mapping as Minimap2_target from './NextflowModules/Minimap2/2.2.4--h5bf99c6_0/Mapping.nf' params(optional: " -y -ax map-ont", genome_fasta: params.target_genome_fasta)
include Mapping as Minimap2_remap from './NextflowModules/Minimap2/2.2.4--h5bf99c6_0/Mapping.nf' params(optional: " -y -ax map-ont", genome_fasta: params.genome_fasta)
include ViewSort as Sambamba_ViewSort_target from './NextflowModules/Sambamba/0.7.0/ViewSort.nf'
include ViewSort as Sambamba_ViewSort_remap from './NextflowModules/Sambamba/0.7.0/ViewSort.nf'


//include HaplotypeCaller_SMN as GATK_HaplotypeCaller_Paraphase from './NextflowModules/GATK/4.2.1.0/HaplotypeCaller.nf' params(genome: params.genome_fasta, compress:true, extention: "_paraphase", optional:"--intervals $params.calling_target_paraphase")
include HaplotypeCaller_SMN as GATK_HaplotypeCaller_Paraphase from './NextflowModules/GATK/4.2.1.0/HaplotypeCaller.nf' params(genome: params.genome_fasta, compress:true, extention: "_paraphase", optional:"--intervals $params.calling_target_paraphase --dont-use-soft-clipped-bases")

include FilterVcfs as GATK_FilterSNV_Paraphase from './NextflowModules/GATK/4.2.1.0/FilterVCFs.nf' params(genome: params.genome_fasta, filter: "SNP")
include Phase as Whatshap_Phase_Target_Paraphase from './NextflowModules/Whatshap/1.7/Phase.nf' params (genome: params.genome_fasta)
include Haplotag as Whatshap_Haplotag_Target_Paraphase from './NextflowModules/Whatshap/1.7/Haplotag.nf' params (genome: params.genome_fasta, extention: "_paraphase")
include Zip_Index as Tabix_Zip_Index_Paraphase from './NextflowModules/Tabix/1.11/Index.nf'
include Index as Sambamba_Index_Target_Paraphase from './NextflowModules/Sambamba/0.7.0/Index.nf'


//include HaplotypeCaller_SMN as GATK_HaplotypeCaller_Region from './NextflowModules/GATK/4.2.1.0/HaplotypeCaller.nf' params(genome: params.genome_fasta, compress:true, extention: "_region", optional:"--intervals $params.calling_target_region")
include HaplotypeCaller_SMN as GATK_HaplotypeCaller_Region from './NextflowModules/GATK/4.2.1.0/HaplotypeCaller.nf' params(genome: params.genome_fasta, compress:true, extention: "_region", optional:"--intervals $params.calling_target_region --dont-use-soft-clipped-bases")

include Phase as Whatshap_Phase_Target_Region from './NextflowModules/Whatshap/1.7/Phase.nf' params (genome: params.genome_fasta)
include Haplotag as Whatshap_Haplotag_Target_Region from './NextflowModules/Whatshap/1.7/Haplotag.nf' params (genome: params.genome_fasta, extention: "_region")
include Zip_Index as Tabix_Zip_Index_Region from './NextflowModules/Tabix/1.11/Index.nf'
include Index as Sambamba_Index_Target_Region from './NextflowModules/Sambamba/0.7.0/Index.nf'


include CollectMultipleMetrics as PICARD_CollectMultipleMetrics from './NextflowModules/Picard/2.22.0/CollectMultipleMetrics.nf' params(genome:"$params.genome_fasta", optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE")
include CollectWgsMetrics as PICARD_CollectWgsMetrics from './NextflowModules/Picard/2.22.0/CollectWgsMetrics.nf' params(genome:"$params.genome_fasta", optional: "MINIMUM_MAPPING_QUALITY=1 MINIMUM_BASE_QUALITY=1 ")
include MultiQC from './NextflowModules/MultiQC/1.9/MultiQC.nf' params(optional: "--config $baseDir/assets/multiqc_config.yaml")

def analysis_id = params.outdir.split('/')[-1]
sample_id = params.sample_id

workflow {

    if( params.start == 'bam'){
        // Get fast5 and mapped bams from input folder
        fast5_files = Channel.fromPath(params.input_path +  "/workspace/fast5_pass/*.fast5").toList()
        bam_files = Channel.fromPath(params.input_path +  "/pass/*.bam").toList()
        summary_file = Channel.fromPath(params.input_path +  "/sequencing_summary.txt").toList()
    }
    else if( params.start == 'rebase' ){
        //Re-basecalling
        ReBasecallingGuppy(params.input_path, sample_id)
        fast5_files = ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files, summary_file -> fast5_files}
        bam_files = ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files, summary_file -> bam_files}
        summary_file = ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files, summary_file -> summary_file}
    }
    else if( params.start == 'bam_remap' ){
        // Get fast5 and mapped bams from input folder
        fast5_files = Channel.fromPath(params.input_path +  "/workspace/fast5_pass/*.fast5").toList()
        bam_files_guppy = Channel.fromPath(params.input_path +  "/pass/*.bam").toList()
        summary_file = Channel.fromPath(params.input_path +  "/sequencing_summary.txt").toList()

        // Extract FASTQ from BAM
        Samtools_Fastq(bam_files_guppy.flatten())

         // Re-map ROI fastq
        Minimap2_remap(Samtools_Fastq.out)

        // Sort SAM to BAM
        Sambamba_ViewSort_remap(Minimap2_remap.out)
 
        bam_files = Sambamba_ViewSort_remap.out.map{sample_id, bam_file, bai_file -> bam_file}
    }
    else{
        error "Invalid alignment mode: ${start}. This should be bam (start from basecalled data), bam_remap (start from bam, but perform remapping with minimap2), or rebase (full re-basecalling)"
    }


    // Add readgroup to BAMs
    Samtools_AddReadgroup(sample_id, bam_files.flatten())

    // Filter for minimum readlength
    Sambamba_Filter_Condition(Samtools_AddReadgroup.out)

    // MergeSort BAMs
    Sambamba_Merge(sample_id, Sambamba_Filter_Condition.out
        .map{bam_file, bai_file -> bam_file}.collect()
    )

    // Identify readpairs
    Duplex_PairsFromSummary(sample_id, summary_file)

    // Identify possible duplex reads from read pairs
    Duplex_FilterPairs(Duplex_PairsFromSummary.out, Sambamba_Merge.out)

    //Filter BAM for duplicate duplex read
    PICARD_FilterSamReads(Sambamba_Merge.out, Duplex_FilterPairs.out)

    //Index BAM file
    Sambamba_Index_Deduplex(PICARD_FilterSamReads.out)

    Bam_file = PICARD_FilterSamReads.out.combine(Sambamba_Index_Deduplex.out
        .map{bam_file, bai_file -> bai_file}
    ).map{bam_file, bai_file -> [sample_id, bam_file, bai_file]}

    if (params.method == "wgs"){
        //Phasing BAM
        LongshotPhase(Bam_file)

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())
    }

    if (params.method == "wgs_repeat"){

        //Phasing BAM
        LongshotPhase(Bam_file)

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())
     
        //SplitPhasedBam
        Sambamba_Split(Sambamba_Index_Longshot.out)

        //Convert BAM to SAM
        Sambamba_ToSam(Sambamba_Index_Longshot.out)
        Sambamba_ToSam_hap1(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap1_bam, hap1_bai]})
        Sambamba_ToSam_hap2(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap2_bam, hap2_bai]})
        Sambamba_ToSam_nohap(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, nohap_bam, nohap_bai]})

        // Get ReadIDs
        Sambamba_GetReadIDs(Sambamba_Index_Longshot.out)
        Sambamba_GetReadIDs_hap1(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap1_bam, hap1_bai]})
        Sambamba_GetReadIDs_hap2(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap2_bam, hap2_bai]})
        Sambamba_GetReadIDs_nohap(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, nohap_bam, nohap_bai]})

        //STRique index
        STRiqueIndex(fast5_files.flatten())

        //Concat all fofn files
        ConcatFofn(STRiqueIndex.out.collect(), sample_id)

        //Repeat calling
        STRiqueCallRepeat(Sambamba_ToSam.out, ConcatFofn.out, fast5_files.collect())
        STRiqueCallRepeat_hap1(Sambamba_ToSam_hap1.out, ConcatFofn.out, fast5_files.collect())
        STRiqueCallRepeat_hap2(Sambamba_ToSam_hap2.out, ConcatFofn.out, fast5_files.collect())
        STRiqueCallRepeat_nohap(Sambamba_ToSam_nohap.out, ConcatFofn.out, fast5_files.collect())

    }
  
    if (params.method == "wgs_roi"){
        // Filter BAM on roi
        Sambamba_Filter_ROI(Bam_file.map{bam_file, bai_file -> bam_file})

        //Phasing roi BAM
        LongshotPhase(Sambamba_Filter_ROI.out)

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())
    }

    if (params.method == "wgs_roi_repeat"){
        // Filter BAM on roi
        Sambamba_Filter_ROI(Bam_file.map{bam_file, bai_file -> bam_file})

        //Phasing roi BAM
        LongshotPhase(Sambamba_Filter_ROI.out)

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())

        //SplitPhasedBam
        Sambamba_Split(Sambamba_Index_Longshot.out)

        //Convert BAM to SAM
        Sambamba_ToSam(Sambamba_Index_Longshot.out)
        Sambamba_ToSam_hap1(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap1_bam, hap1_bai]})
        Sambamba_ToSam_hap2(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap2_bam, hap2_bai]})
        Sambamba_ToSam_nohap(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, nohap_bam, nohap_bai]})

        // Get ReadIDs
        Sambamba_GetReadIDs(Sambamba_Index_Longshot.out)
        Sambamba_GetReadIDs_hap1(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap1_bam, hap1_bai]})
        Sambamba_GetReadIDs_hap2(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap2_bam, hap2_bai]})
        Sambamba_GetReadIDs_nohap(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, nohap_bam, nohap_bai]})

        //STRique index
        STRiqueIndex(fast5_files.flatten())

        //Concat all fofn files
        ConcatFofn(STRiqueIndex.out.collect(), sample_id)

        //Repeat calling
        STRiqueCallRepeat(Sambamba_ToSam.out, ConcatFofn.out, fast5_files.collect())
        STRiqueCallRepeat_hap1(Sambamba_ToSam_hap1.out, ConcatFofn.out, fast5_files.collect())
        STRiqueCallRepeat_hap2(Sambamba_ToSam_hap2.out, ConcatFofn.out, fast5_files.collect())
        STRiqueCallRepeat_nohap(Sambamba_ToSam_nohap.out, ConcatFofn.out, fast5_files.collect())

    }

    if (params.method == "wgs_splitcas9_repeat"){
        // BAM split based on Cas9 sites
        SplitBAM(Bam_file)

        ParseSampleIDs = Channel.fromPath( params.splitfile )
            .splitCsv( sep: '\t' )
            .map{sample_id, chromosome, start, stop -> [sample_id]}
            .unique()
        
        //Phasing BAMs
        LongshotPhase(SplitBAM.out.transpose()
            .map { tuple(it) }
            .map{sample_id, bam_file, bai_file -> [bam_file.simpleName.toString().split("_")[0], bam_file, bai_file]}
            .join(ParseSampleIDs)
        )

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())

        //SplitPhasedBam
        Sambamba_Split(Sambamba_Index_Longshot.out)

        //Convert BAM to SAM
        Sambamba_ToSam(Sambamba_Index_Longshot.out)
        Sambamba_ToSam_hap1(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap1_bam, hap1_bai]})
        Sambamba_ToSam_hap2(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap2_bam, hap2_bai]})
        Sambamba_ToSam_nohap(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, nohap_bam, nohap_bai]})

        // Get ReadIDs
        Sambamba_GetReadIDs(Sambamba_Index_Longshot.out)
        Sambamba_GetReadIDs_hap1(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap1_bam, hap1_bai]})
        Sambamba_GetReadIDs_hap2(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap2_bam, hap2_bai]})
        Sambamba_GetReadIDs_nohap(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, nohap_bam, nohap_bai]})

        //STRique index
        STRiqueIndex(fast5_files.flatten())

        //Concat all fofn files
        ConcatFofn(STRiqueIndex.out.collect(), sample_id)
    
        //Repeat calling
        STRiqueCallRepeat(Sambamba_ToSam.out, ConcatFofn.out, fast5_files.collect())
        STRiqueCallRepeat_hap1(Sambamba_ToSam_hap1.out, ConcatFofn.out, fast5_files.collect())
        STRiqueCallRepeat_hap2(Sambamba_ToSam_hap2.out, ConcatFofn.out, fast5_files.collect())
        STRiqueCallRepeat_nohap(Sambamba_ToSam_nohap.out, ConcatFofn.out, fast5_files.collect())
    }

    if (params.method == "targeted_splitcas9"){

        // Make FASTQ files of ROI including tags
        Samtools_Fastq(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [bam_file, bai_file]}) 
 
        // Re-map ROI fastq
        Minimap2_target(Samtools_Fastq.out)

        // Sort SAM to BAM
        Sambamba_ViewSort_target(Minimap2_target.out)

        // BAM split on cas9 sites
        SplitBAM(Sambamba_ViewSort_target.out)
        //SplitBAM(Bam_file)

        ParseSampleIDs = Channel.fromPath( params.splitfile )
            .splitCsv( sep: '\t' )
            .map{sample_id, chromosome, start, stop -> [sample_id]}
            .unique()

        //Phasing BAM
        LongshotPhase(SplitBAM.out.transpose()
            .map { tuple(it) }
            .map{sample_id, bam_file, bai_file -> [bam_file.simpleName.toString().split("_")[0], bam_file, bai_file]}
            .join(ParseSampleIDs)
        )

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())

    }

    
    if (params.method == "targeted_SMA_splitcas9"){

        // Make FASTQ files of ROI including tags
        Samtools_Fastq(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [bam_file, bai_file]})

        // Re-map ROI fastq
        Minimap2_target(Samtools_Fastq.out)

        // Sort SAM to BAM
        Sambamba_ViewSort_target(Minimap2_target.out)

        // BAM split based on Cas9 sites
        //SplitBAM(Sambamba_Merge.out)
        SplitBAM(Bam_file)

        // Variant calling
        ParsePloidy = Channel.fromPath( params.splitfile )
            .splitCsv( sep: '\t' )
            .map{sample_id, chromosome, start, stop, ploidy -> [sample_id, ploidy]}
            .unique()
 
        GATK_HaplotypeCaller_Paraphase(SplitBAM.out.transpose()
            .map { tuple(it) }
            .map{sample_id, bam_file, bai_file -> [bam_file.simpleName.toString().split("_")[0], bam_file, bai_file]}
            .join(ParsePloidy)
        )

        // Filter SNV only
        GATK_FilterSNV(GATK_HaplotypeCaller_Paraphase.out)

        // Whatshapp polyphase 
        Whatshap_Phase_Target(GATK_FilterSNV.out)
     
        // bgzip and index VCF
        Tabix_Zip_Index(Whatshap_Phase_Target.out)
     
        // Whatshapp haplotag
        Whatshap_Haplotag_Target(
            GATK_HaplotypeCaller_Paraphase.out
            .map{sample_id, bam_file, bai_file, vcf_file, vcf_index, ploidy -> [sample_id, bam_file, bai_file]}
            .join(Tabix_Zip_Index.out)
        )

        // Index BAM file and publish
        Sambamba_Index_Target(Whatshap_Haplotag_Target.out)
    }

    if (params.method == "targeted_SMA_adaptive"){
        GATK_HaplotypeCaller_Paraphase(Bam_file
             .map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file, params.ploidy]})

        // Filter SNV only Paraphase variants
        GATK_FilterSNV_Paraphase(GATK_HaplotypeCaller_Paraphase.out)

        // Whatshapp polyphase Paraphase variants
        Whatshap_Phase_Target_Paraphase(GATK_FilterSNV_Paraphase.out)

        // bgzip and index VCF Paraphase variants
        Tabix_Zip_Index_Paraphase(Whatshap_Phase_Target_Paraphase.out)

        // Whatshapp haplotag Paraphase variants
        Whatshap_Haplotag_Target_Paraphase(
            GATK_HaplotypeCaller_Paraphase.out
            .map{sample_id, bam_file, bai_file, vcf_file, vcf_index, ploidy -> [sample_id, bam_file, bai_file]}
            .join(Tabix_Zip_Index_Paraphase.out)
        )

        // Index BAM file and publish Paraphase variants
        Sambamba_Index_Target_Paraphase(Whatshap_Haplotag_Target_Paraphase.out)


        // Variant calling on used defined region
        GATK_HaplotypeCaller_Region(Bam_file
             .map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file, params.ploidy]})

        // Whatshapp polyphase region
        Whatshap_Phase_Target_Region(GATK_HaplotypeCaller_Region.out
            .map{sample_id, bam_file, bai_file, vcf, vcf_index, ploidy -> [sample_id, bam_file, bai_file, vcf, ploidy]}
        )

        // bgzip and index VCF
        Tabix_Zip_Index_Region(Whatshap_Phase_Target_Region.out)

        // Whatshapp haplotag
        Whatshap_Haplotag_Target_Region(
            GATK_HaplotypeCaller_Region.out
            .map{sample_id, bam_file, bai_file, vcf_file, vcf_index, ploidy -> [sample_id, bam_file, bai_file]}
            .join(Tabix_Zip_Index_Region.out)
        )

        // Index BAM file and publish region variants
        Sambamba_Index_Target_Region(Whatshap_Haplotag_Target_Region.out)
    }


    // QC stats
    PICARD_CollectMultipleMetrics(Bam_file)
    PICARD_CollectWgsMetrics(Bam_file)
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

// Custom processes
process ReBasecallingGuppy{
    // Custom process to re-basecall data with guppy
    tag {"ReBasecallingGuppy ${sample_id}"}
    label 'ReBasecallingGuppy'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        val(input_path)
        val(sample_id)

    output:
        tuple(path("pass/*.fastq.gz"), path("workspace/fast5_pass/*.fast5"), path ("*"), path("pass/*.bam"), path("sequencing_summary.txt"))
        //tuple(path("pass/*.fastq.gz"), path("workspace/*.fast5"), path ("*"), path("pass/*.bam"))

    script:
        // adding --index will also give .bai in output. Not implemented yet.     
        """
        $params.guppy_basecaller_path -x "cuda:0" -c $params.guppy_path/data/$params.guppy_basecaller_config \
        -i $input_path -s ./ $params.guppy_basecaller_params \
        --bam_out --fast5_out --align_ref $params.genome_mapping_index 
        """
}

process SplitBAM{
    // Custom process to split BAM based on Cas9 start sites
    tag {"SplitBAM ${sample_id}"}
    label 'SplitBAM'
    shell = ['/bin/bash', '-eo', 'pipefail']
    //cache = false

    input:
        tuple(sample_id, path(bam_file), path(bai_file))

    output:
        tuple(sample_id, "*split.bam", "*split.bam.bai")

    script:
        """
        source /hpc/diaggen/users/Martin/Research_projects_Martin/SMA_project/venv/bin/activate
        python /hpc/diaggen/users/Martin/Research_projects_Martin/SMA_project/split_bam_start_site.py $bam_file $params.splitfile --flanks 20
        """
}


process ConcatFofn{
    // Custom process to concat all fofn files
    tag {"Concat Fofn ${sample_id}"}
    label 'Concat_Fofn'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        path(fofn_files)
        val(sample_id)

    output:
        path("${sample_id}.fofn")

    script:
        fofn  = fofn_files.join(' ')
        """
        cat ${fofn} > ${sample_id}.fofn
        """
}


process VersionLog {
    // Custom process to log repository versions
    tag {"VersionLog ${analysis_id}"}
    label 'VersionLog'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a new version log when restarting the workflow.

    output:
        path('repository_version.log')

    script:
        """
        echo 'DxNextflowONT' > repository_version.log
        git --git-dir=${workflow.projectDir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        """
}

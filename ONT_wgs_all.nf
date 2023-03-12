#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Utils modules
include ExportParams as Workflow_ExportParams from './NextflowModules/Utils/workflow.nf'

include Index as Sambamba_Index from './NextflowModules/Sambamba/0.7.0/Index.nf'
include Merge as Sambamba_Merge from './NextflowModules/Sambamba/0.7.0/Merge.nf'
include Filter as Sambamba_Filter from './NextflowModules/Sambamba/0.7.0/Filter.nf'
include Split as Sambamba_Split from './NextflowModules/Sambamba/0.7.0/Split.nf'
include Index as Sambamba_Index_Longshot from './NextflowModules/Sambamba/0.7.0/Index.nf'
include LongshotPhase from './NextflowModules/Longshot/0.4.1/Phase.nf'

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
include Mapping as Minimap2_mapping from "./NextflowModules/Minimap2/2.2.4--h5bf99c6_0/Mapping.nf" params(optional: " -y -ax map-ont", genome_fasta: params.target_genome_fasta)
include ViewSort as Sambamba_ViewSort from "./NextflowModules/Sambamba/0.7.0/ViewSort.nf"
include HaplotypeCaller_SMN as GATK_HaplotypeCaller_SMN from "./NextflowModules/GATK/4.2.1.0/HaplotypeCaller.nf" params(genome: params.target_genome_fasta, compress:true, optional:"")
include FilterVcfs as GATK_FilterSNV from "./NextflowModules/GATK/4.2.1.0/FilterVCFs.nf" params(genome: params.target_genome_fasta, filter: "SNP")
include Phase as Whatshap_Phase_Target from "./NextflowModules/Whatshap/1.7/Phase.nf" params (genome: params.target_genome_fasta)

include CollectMultipleMetrics as PICARD_CollectMultipleMetrics from './NextflowModules/Picard/2.22.0/CollectMultipleMetrics.nf' params(genome:"$params.genome", optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE")
include CollectWgsMetrics as PICARD_CollectWgsMetrics from './NextflowModules/Picard/2.22.0/CollectWgsMetrics.nf' params(genome:"$params.genome", optional: "MINIMUM_MAPPING_QUALITY=1 MINIMUM_BASE_QUALITY=1 ")
include MultiQC from './NextflowModules/MultiQC/1.9/MultiQC.nf' params(optional: "--config $baseDir/assets/multiqc_config.yaml")

def analysis_id = params.outdir.split('/')[-1]
sample_id = params.sample_id

workflow {

    //Re-basecalling
    ReBasecallingGuppy(params.fast5_path, sample_id)

    //BAMindex
    Sambamba_Index(sample_id, ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files -> bam_files}.flatten())

    // MergeSort BAMs
    Sambamba_Merge(Sambamba_Index.out.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}.groupTuple())

    if (params.method == "wgs"){

        //Phasing BAM
        LongshotPhase(Sambamba_Merge.out)

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())
    }
  
    if (params.method == "wgs_roi"){
 
        // Filter BAM on roi
        Sambamba_Filter(Sambamba_Merge.out)

        //Phasing roi BAM
        LongshotPhase(Sambamba_Filter.out)

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())
    }

    if (params.method == "wgs_roi_repeat"){
        // Filter BAM on roi
        Sambamba_Filter(Sambamba_Merge.out)

        //Phasing roi BAM
        LongshotPhase(Sambamba_Filter.out)

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())

        //SplitPhasedBam
        Sambamba_Split(LongshotPhase.out)

        //Convert BAM to SAM
        Sambamba_ToSam(Sambamba_Filter.out)
        Sambamba_ToSam_hap1(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap1_bam, hap1_bai]})
        Sambamba_ToSam_hap2(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap2_bam, hap2_bai]})
        Sambamba_ToSam_nohap(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, nohap_bam, nohap_bai]})

        // Get ReadIDs
        Sambamba_GetReadIDs(Sambamba_Filter.out)
        Sambamba_GetReadIDs_hap1(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap1_bam, hap1_bai]})
        Sambamba_GetReadIDs_hap2(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, hap2_bam, hap2_bai]})
        Sambamba_GetReadIDs_nohap(Sambamba_Split.out.map{sample_id, hap1_bam, hap1_bai, hap2_bam, hap2_bai, nohap_bam, nohap_bai -> [sample_id, nohap_bam, nohap_bai]})

        //STRique index
        STRiqueIndex(ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files -> fast5_files}.flatten())

        //Concat all fofn files
        ConcatFofn(STRiqueIndex.out.collect(), sample_id)

    }

    if (params.method == "wgs_splitcas9"){

        //Phasing BAM
        LongshotPhase(Sambamba_Merge.out)

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())

        // BAM split based on Cas9 sites
        SplitBAM(Sambamba_Merge.out)
        //SplitBAM.out.transpose().map { tuple(it) }.view()
    }

    if (params.method == "targeted"){

        // Make FASTQ files of ROI including tags
        Samtools_Fastq(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [bam_file, bai_file]})

        // Re-map ROI fastq
        Minimap2_mapping(Samtools_Fastq.out)

        // Sort SAM to BAM
        Sambamba_ViewSort(Minimap2_mapping.out)  

    }


    if (params.method == "targeted_splitcas9"){

        // Make FASTQ files of ROI including tags
        Samtools_Fastq(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [bam_file, bai_file]}) 
 
        // Re-map ROI fastq
        Minimap2_mapping(Samtools_Fastq.out)

        // Sort SAM to BAM
        Sambamba_ViewSort(Minimap2_mapping.out)

        // BAM split on cas9 sites
        SplitBAM(Sambamba_Merge.out)

        // Make FASTQ files of ROI including tags
        Samtools_Fastq(SplitBAM.out.transpose().map { tuple(it) })

    }

    if (params.method == "targeted_SMA_splitcas9"){

        // Make FASTQ files of ROI including tags
        Samtools_Fastq(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [bam_file, bai_file]})

        // Re-map ROI fastq
        Minimap2_mapping(Samtools_Fastq.out)

        // Sort SAM to BAM
        Sambamba_ViewSort(Minimap2_mapping.out)

        // BAM split based on Cas9 sites
        SplitBAM(Sambamba_ViewSort.out)

        // Variant calling
        ParsePloidy = Channel.fromPath( params.splitfile )
            .splitCsv( sep: '\t' )
            .map{sample_id, chromosome, start, stop, ploidy -> [sample_id, ploidy]}
  
        GATK_HaplotypeCaller_SMN(SplitBAM.out.transpose()
            .map { tuple(it) }
            .map{sample_id, bam_file, bai_file -> [bam_file.simpleName.toString().split("_")[0], bam_file, bai_file]}
            .join(ParsePloidy)
        )

        //GATK_HaplotypeCaller_SMN.out.view()

        // Filter SNV only
        GATK_FilterSNV(GATK_HaplotypeCaller_SMN.out)
        GATK_FilterSNV.out.view()

        // Whatshapp polyphase 
        Whatshap_Phase_Target(GATK_FilterSNV.out)
        // Whatshapp haplotag

        // Haplotag BAM based on VC
        // Haplotag(HaplotypeCaller.out)

    }


    if (params.method == "targeted_SMA_adaptive"){

        // Make FASTQ files of ROI including tags
        Samtools_Fastq(Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [bam_file, bai_file]})

        // Re-map ROI fastq
        Minimap2_mapping(Samtools_Fastq.out)
 
        // Sort SAM to BAM
        Sambamba_ViewSort(Minimap2_mapping.out)

        // Variant calling
        //GATK_HaplotypeCaller_SMN(Sambamba_ViewSort.out)       

        // Haplotag BAM based on VC
        // Haplotag(HaplotypeCaller.out)


    }


    // QC stats
    //PICARD_CollectMultipleMetrics(Sambamba_Merge.out)
    //PICARD_CollectWgsMetrics(Sambamba_Merge.out)
    //MultiQC(analysis_id, Channel.empty().mix(
    //    PICARD_CollectMultipleMetrics.out,
    //    PICARD_CollectWgsMetrics.out
    //).collect())

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
        val(fast5_path)
        val(sample_id)

    output:
        tuple(path("pass/*.fastq.gz"), path("workspace/fast5_pass/*.fast5"), path ("*"), path("pass/*.bam"))
        //tuple(path("pass/*.fastq.gz"), path("workspace/*.fast5"), path ("*"), path("pass/*.bam"))

    script:
        """
        $params.guppy_basecaller_path -x "cuda:0" -c $params.guppy_path/data/$params.guppy_basecaller_config \
        --num_callers ${task.cpus} -i $fast5_path -s ./ $params.guppy_basecaller_params \
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
        tuple(sample_id, rg_id, path(bam_file), path(bai_file))

    output:
        tuple(sample_id, "*split.bam", "*split.bam.bai")

    script:
        """
        source /hpc/diaggen/users/Martin/Research_projects_Martin/SMA_project/venv/bin/activate
        python /hpc/diaggen/users/Martin/Research_projects_Martin/SMA_project/split_bam_start_site.py $bam_file $params.splitfile --flanks 20
        """
}


process SaveInputFile {
    tag {"SaveInputFile ${analysis_id}"}
    label 'SaveInputFile'
    shell = ['/bin/bash', '-euo', 'pipefail']
    cache = false  //Disable cache to force a new files to be copied.
  
    input:
       path(roi)
 
    output:
       path(roi)

    script:
        """
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

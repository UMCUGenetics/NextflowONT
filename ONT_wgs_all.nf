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

include Fastq as Samtools_Fastq from './NextflowModules/Samtools/1.15/Fastq.nf' params(tags: " -T RG,Mm,Ml ")
include Mapping as Minimap2_mapping from "./NextflowModules/Minimap2/2.2.4--h5bf99c6_0/Mapping.nf" params(optional: " -y -ax map-ont", genome_fasta: params.target_genome_fasta)
include ViewSort as Sambamba_ViewSort from "./NextflowModules/Sambamba/0.7.0/ViewSort.nf"

include CollectMultipleMetrics as PICARD_CollectMultipleMetrics from './NextflowModules/Picard/2.22.0/CollectMultipleMetrics.nf' params(genome:"$params.genome", optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE")
include CollectWgsMetrics as PICARD_CollectWgsMetrics from './NextflowModules/Picard/2.22.0/CollectWgsMetrics.nf' params(genome:"$params.genome", optional: "MINIMUM_MAPPING_QUALITY=1 MINIMUM_BASE_QUALITY=1 ")
include MultiQC from './NextflowModules/MultiQC/1.9/MultiQC.nf' params(optional: "--config $baseDir/assets/multiqc_config.yaml")

def analysis_id = params.outdir.split('/')[-1]
sample_id = params.sample_id

workflow {

    if (params.method == "wgs"){
        //Re-basecalling
        ReBasecallingGuppy(params.fast5_path, sample_id)

        //BAMindex
        Sambamba_Index(sample_id, ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files -> bam_files}.flatten())

        // MergeSort BAMs
        Sambamba_Merge(Sambamba_Index.out.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}.groupTuple())

        //Phasing BAM
        LongshotPhase(Sambamba_Merge.out)

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())

        // QC stats
        PICARD_CollectMultipleMetrics(Sambamba_Merge.out)
        PICARD_CollectWgsMetrics(Sambamba_Merge.out)

        MultiQC(analysis_id, Channel.empty().mix(
            PICARD_CollectMultipleMetrics.out,
            PICARD_CollectWgsMetrics.out
        ).collect())
    }

    if (params.method == "targeted"){
        //Re-basecalling
        ReBasecallingGuppy(params.fast5_path, sample_id)

        //BAMindex
        Sambamba_Index(sample_id, ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files -> bam_files}.flatten())

        // MergeSort BAMs
        Sambamba_Merge(Sambamba_Index.out.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}.groupTuple())

        //Phasing BAM
        LongshotPhase(Sambamba_Merge.out)

        // BAMIndex
        Sambamba_Index_Longshot(sample_id, LongshotPhase.out.map{sample_id, bam_file, vcf_file -> bam_file}.flatten())

        // Make FASTQ files of ROI including tags
        Samtools_Fastq(Sambamba_Merge.out)

        // Re-map ROI fastq
        Minimap2_mapping(Samtools_Fastq.out)

        // Sort SAM to BAM
        Sambamba_ViewSort(Minimap2_mapping.out)  

        // QC stats
        PICARD_CollectMultipleMetrics(Sambamba_Merge.out)
        PICARD_CollectWgsMetrics(Sambamba_Merge.out)

        MultiQC(analysis_id, Channel.empty().mix(
            PICARD_CollectMultipleMetrics.out,
            PICARD_CollectWgsMetrics.out
        ).collect())


    }

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
        //tuple(path("pass/*.fastq.gz"), path("workspace/fast5_pass/*.fast5"), path ("*"), path("pass/*.bam"))
        tuple(path("pass/*.fastq.gz"), path("workspace/*.fast5"), path ("*"), path("pass/*.bam"))

    script:
        """
        $params.guppy_basecaller_path -x "cuda:0" -c $params.guppy_path/data/$params.guppy_basecaller_config \
        --num_callers ${task.cpus} -i $fast5_path -s ./ $params.guppy_basecaller_params \
        --bam_out --fast5_out --align_ref $params.genome_mapping_index 
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

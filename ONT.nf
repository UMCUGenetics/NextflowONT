#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Utils modules
include ExportParams as Workflow_ExportParams from './NextflowModules/Utils/workflow.nf'

include Index as Sambamba_Index from './NextflowModules/Sambamba/0.7.0/Index.nf'
include Merge as Sambamba_Merge from './NextflowModules/Sambamba/0.7.0/Merge.nf'
include Filter as Sambamba_Filter from './NextflowModules/Sambamba/0.7.0/Filter.nf'
include Split as Sambamba_Split from './NextflowModules/Sambamba/0.7.0/Split.nf'

// ## Combine to one module. Include optional parameters to or emit diffent hap groups?
include ToSAM as Sambamba_ToSam from './NextflowModules/Sambamba/0.7.0/ToSAM.nf'
include ToSAM as Sambamba_ToSam_hap1 from './NextflowModules/Sambamba/0.7.0/ToSAM.nf'
include ToSAM as Sambamba_ToSam_hap2 from './NextflowModules/Sambamba/0.7.0/ToSAM.nf'
include ToSAM as Sambamba_ToSam_nohap from './NextflowModules/Sambamba/0.7.0/ToSAM.nf'

// ## Combine to one module. Include optional parameters to or emit diffent hap groups?
include GetReadIDs as Sambamba_GetReadIDs from './NextflowModules/Sambamba/0.7.0/GetReads.nf'
include GetReadIDs as Sambamba_GetReadIDs_hap1 from './NextflowModules/Sambamba/0.7.0/GetReads.nf'
include GetReadIDs as Sambamba_GetReadIDs_hap2 from './NextflowModules/Sambamba/0.7.0/GetReads.nf'
include GetReadIDs as Sambamba_GetReadIDs_nohap from './NextflowModules/Sambamba/0.7.0/GetReads.nf'

include LongshotPhase from './NextflowModules/Longshot/0.4.1/Phase.nf'
include Index as STRiqueIndex from './NextflowModules/STRique/0.4.2/Index.nf'

// ## Combine to one module. Include optional parameters to or emit diffent hap groups?
include CallRepeat as STRiqueCallRepeat from './NextflowModules/STRique/0.4.2/CallRepeats.nf'
include CallRepeat as STRiqueCallRepeat_hap1 from './NextflowModules/STRique/0.4.2/CallRepeats.nf'
include CallRepeat as STRiqueCallRepeat_hap2 from './NextflowModules/STRique/0.4.2/CallRepeats.nf'
include CallRepeat as STRiqueCallRepeat_nohap from './NextflowModules/STRique/0.4.2/CallRepeats.nf'

def analysis_id = params.outdir.split('/')[-1]
sample_id = params.sample_id

workflow {
    //Re-basecalling
    ReBasecallingGuppy(params.fast5_path, sample_id)

    //BAMindex
    Sambamba_Index(sample_id, ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files -> bam_files}.flatten())

    // MergeSort BAMs
    Sambamba_Merge(Sambamba_Index.out.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}.groupTuple())
 
    // Filter BAM on roi
    Sambamba_Filter(Sambamba_Merge.out)

    //Phasing BAM
    LongshotPhase(Sambamba_Filter.out) 

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
    
    //Repeat calling
    STRiqueCallRepeat(Sambamba_ToSam.out, ConcatFofn.out, ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files -> fast5_files}.collect())
    STRiqueCallRepeat_hap1(Sambamba_ToSam_hap1.out, ConcatFofn.out, ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files -> fast5_files}.collect())
    STRiqueCallRepeat_hap2(Sambamba_ToSam_hap2.out, ConcatFofn.out, ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files -> fast5_files}.collect())
    STRiqueCallRepeat_nohap(Sambamba_ToSam_nohap.out, ConcatFofn.out, ReBasecallingGuppy.out.map{fastq_files, fast5_files, all_files, bam_files -> fast5_files}.collect())

    //Save Input files
    SaveInputFile(params.roi)

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
        //tuple(path("pass/*.fastq.gz"), path("workspace/*.fast5"), path ("*"))

    script:
        """
        $params.guppy_basecaller_path -x "cuda:0" -c $params.guppy_path/data/$params.guppy_basecaller_config \
        --num_callers ${task.cpus} -i $fast5_path -s ./ $params.guppy_basecaller_params \
        --bam_out --fast5_out --align_ref $params.genome_mapping_index 
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
        //println fofn_files[0]
        """
        cat ${fofn} > ${sample_id}.fofn
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

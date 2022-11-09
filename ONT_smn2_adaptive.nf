#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Utils modules
include ExportParams as Workflow_ExportParams from './NextflowModules/Utils/workflow.nf'

include Index as Sambamba_Index from './NextflowModules/Sambamba/0.7.0/Index.nf'
include Merge as Sambamba_Merge from './NextflowModules/Sambamba/0.7.0/Merge.nf'
include Filter as Sambamba_Filter from './NextflowModules/Sambamba/0.7.0/Filter.nf'

def analysis_id = params.outdir.split('/')[-1]
sample_id = params.sample_id

workflow {
    //Re-basecalling
    ReBasecallingGuppy(params.fast5_path, sample_id)

    //BAMindex
    Sambamba_Index(sample_id, ReBasecallingGuppy.out.map{fastq_files, all_files, bam_files -> bam_files}.flatten())

    // MergeSort BAMs
    Sambamba_Merge(Sambamba_Index.out.map{sample_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}.groupTuple())
 
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
        tuple(path("pass/*.fastq.gz"), path ("*"), path("pass/*.bam"))

    script:
        """
        $params.guppy_basecaller_path -x "cuda:0" -c $params.guppy_path/data/$params.guppy_basecaller_config \
        --num_callers ${task.cpus} -i $fast5_path -s ./ $params.guppy_basecaller_params \
        --bam_out --fast5_out --align_ref $params.genome_mapping_index
        """
}

process Split_BAM{
    // Custom process split BAM on (Cas9) start and stop sites
    tag {"SplitBam ${sample_id}"}
    label 'SplitBam'
    shell = ['/bin/bash', '-eo', 'pipefail']
 
    input:
        tuple(val(sample_id), path(bam_file), path(bai_files))
        
    output:
        tuple(path("*.bam"), path("*bai"))

    script:
        """
        source ${params.split_path}/venv/bin/activate
        python ${params.split_path}/split_bam_start_site.py $bam_file $params.split_file --flank $params.split_flanks ${params.split_optional}
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

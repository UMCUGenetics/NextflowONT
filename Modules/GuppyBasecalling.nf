process ReBasecallingGuppy{
    // Custom process to re-basecall data with guppy
    tag {"ReBasecallingGuppy ${sample_id}"}
    label 'ReBasecallingGuppy'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        val(input_path)
        val(sample_id)
        val(fast5)

    output:
        tuple(
            path("pass/*.fastq.gz"),
            path("workspace/fast5_pass/*.fast5"),
            path ("*"),
            path("pass/*.bam"),
            path("sequencing_summary.txt")
        )

    script:
        // adding --index will also give .bai in output. Not implemented yet.
        """
        $params.guppy_basecaller_path -x "cuda:0" -c $params.guppy_path/data/$params.guppy_basecaller_config \
        -i $input_path -s ./ $params.guppy_basecaller_params \
        --bam_out --fast5_out --align_ref $params.genome_mapping_index
        """
}


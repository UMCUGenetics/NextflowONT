process GetPhaseSet {
    tag { "Get_PhaseSet ${vcf_file.baseName}" }
    label 'Get_PhaseSet'
    container = 'quay.io/biocontainers/vcfpy:0.13.6--pyhdfd78af_0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(val(analysis_id), path(vcf_file), path(vcf_index_file))

    output:
        tuple(stdout)

    script:
        """
        $baseDir/assets/get_phaseset.py \
            $vcf_file \
            $params.psv_region | tr -d '\n'
        """
}

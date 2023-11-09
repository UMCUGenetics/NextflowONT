#!/bin/bash
set -euo pipefail

workflow_path='/hpc/diaggen/software/development/NextflowONT'

# Set input and output dirs
input_path=`realpath -e $1`
output=$2
email=$3
sample_id=$4
start_run=$5
method=${6-default}


if [ $method == "wgs" ]; then
    echo " #### Running method wgs  ####"
    optional_params=( "${@:7}" )
fi

if [ $method == "wgs_repeat"  ]; then
    echo " #### Running method wgs_repeat  ####"
    strique_config='--strique_config '$7
    optional_params=( "${@:8}" )
fi

if [ $method == "wgs_roi_repeat" ]; then
    echo " #### Running method wgs + roi phasing + repeat calling ####"
    roi='--roi '$7
    strique_config='--strique_config '$8
    optional_params=( "${@:9}" )
fi

if [ $method == "wgs_splitcas9_repeat" ]; then
    echo " #### Running method wgs + split cas9 ####"
    splitfile='--splitfile '$7
    strique_config='--strique_config '$8
    optional_params=( "${@:9}" )
fi

if [ $method == "SMA_adaptive" ]; then
    echo " #### Running method targeted SMA specific + adaptive sequencing ####"
    ploidy='--ploidy '$7
    optional_params=( "${@:8}" )
fi


mkdir -p $output && cd $output
mkdir -p log

if ! { [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; }; then
touch workflow.running

sbatch <<EOT
#!/bin/bash
#SBATCH --time=144:00:00
#SBATCH --nodes=1
#SBATCH --mem 10G
#SBATCH --gres=tmpspace:16G
#SBATCH --job-name Nextflow_ONT
#SBATCH -o log/slurm_nextflow_ont.%j.out
#SBATCH -e log/slurm_nextflow_ont.%j.err
#SBATCH --mail-user $email
#SBATCH --mail-type FAIL
#SBATCH --export=NONE
#SBATCH --account=diaggen

/hpc/diaggen/software/tools/nextflow run $workflow_path/ONT_wgs_all.nf \
-c $workflow_path/ONT_wgs_all.config \
--input_path $input_path \
--outdir $output \
--email $email \
--sample_id $sample_id \
--start $start_run \
--method $method \
${roi:-""} \
${strique_config:-""} \
${splitfile:-""} \
${ploidy:-""} \
${optional_params[@]:-""} \
-profile slurm \
-resume -ansi-log false

if [ \$? -eq 0 ]; then
    echo "Nextflow done."

    #  echo "Zip work directory"
    find work -type f | egrep "\.(command|exitcode)" | zip -@ -q work.zip

    #  echo "Remove work directory"
    rm -r work

    #  echo "Creating md5sum"
    find -type f -not -iname 'md5sum.txt' -exec md5sum {} \; > md5sum.txt

    echo "WES workflow completed successfully."
    rm workflow.running
    touch workflow.done

    exit 0
else
    echo "Nextflow failed"
    rm workflow.running
    touch workflow.failed
    exit 1
fi
EOT
else
echo "Workflow job not submitted, please check $output for 'workflow.status' files."
fi

#!/bin/bash
set -euo pipefail

workflow_path='/hpc/diaggen/software/development/NextflowONT'

# Set input and output dirs
input_fast5=`realpath -e $1`
output=$2
email=$3
sample_id=$4
method=${5-default}


if [ $method == "wgs" ]; then
    echo " #### Running method wgs  ####"
fi


if [ $method == "wgs_roi" ]; then
    echo " #### Running method wgs + roi phasing  ####"
    roi='--roi '$6
fi


if [ $method == "wgs_roi_repeat" ]; then
    echo " #### Running method wgs + roi phasing + repeat calling ####"
    roi='--roi '$6
    strique_config='--strique_config '$7
fi


if [ $method == "wgs_splitcas9" ]; then
    echo " #### Running method wgs + split cas9 ####"
    splitfile='--splitfile '$7
fi


if [ $method == "targeted" ]; then
    echo " #### Running method wgs + targeted mapping  ####"
    roi='--roi '$6
fi


if [ $method == "targeted_splitcas9" ]; then
    echo " #### Running method targeted + split cas9  ####"
    roi='--roi '$6
    splitfile='--splitfile '$7
fi


if [ $method == "targeted_SMA_splitcas9" ]; then
    echo " #### Running method targeted SMA specific + split cas9  ####"
    roi='--roi '$6
    splitfile='--splitfile '$7
fi


if [ $method == "targeted_SMA_adaptive" ]; then
    echo " #### Running method targeted SMA specific + adaptive sequencing ####"
    roi='--roi '$6
fi


mkdir -p $output && cd $output
mkdir -p log

if ! { [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; }; then
touch workflow.running

sbatch <<EOT
#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --gres=tmpspace:10G
#SBATCH --job-name Nextflow_ONT
#SBATCH -o log/slurm_nextflow_ont.%j.out
#SBATCH -e log/slurm_nextflow_ont.%j.err
#SBATCH --mail-user $email
#SBATCH --mail-type FAIL
#SBATCH --export=NONE
#SBATCH --account=diaggen

module load Java/1.8.0_60

/hpc/diaggen/software/tools/nextflow run $workflow_path/ONT_wgs_all.nf \
-c $workflow_path/ONT_wgs_all.config \
--fast5_path $input_fast5 \
--outdir $output \
--email $email \
--sample_id $sample_id \
--method $method \
${roi:-""} \
${strique_config:-""} \
${splitfile:-""} \
-profile slurm \
-resume -ansi-log false

if [ \$? -eq 0 ]; then
    echo "Nextflow done."

    #  echo "Zip work directory"
    #  find work -type f | egrep "\.(command|exitcode)" | zip -@ -q work.zip

    #  echo "Remove work directory"
    #  rm -r work

    #  echo "Creating md5sum"
    #  find -type f -not -iname 'md5sum.txt' -exec md5sum {} \; > md5sum.txt

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

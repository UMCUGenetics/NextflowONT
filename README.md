# NextflowONT
Nextflow Oxford Nanopore Technologies workflow

#### Get Nextflow Modules
```bash
git submodule update --init --recursive
```

#### Install Nextflow
```bash
mkdir tools && cd tools
curl -s https://get.nextflow.io | bash
```

#### Running ONT workflow
```bash
nextflow run ONT.nf -c ONT.config --fastq_path <fastq_dir_path> --outdir <output_dir_path> --start <bam|rebase> --method <method> --email <email> [-profile slurm|mac]
```
\
\
<methods>
rebase = include re-basecalling\
bam = start from Guppy folder including bam files.\
bam_remap = start from Guppy folder including bam files, but perform remapping to genome in config.
\

| Method | Description | Optional parameters needed|
| --- | :--- | :--- |
|wgs|perform whole genome mapping + longshot phasing|none|
|wgs_repeat|perform whole genome mapping + longshot phasing + repeat calling STRique|--strique_config \<strique_config\>|
wgs_roi|perform whole genome mapping + slice on ROI + longshot phasing|--roi \<roi\>|
wgs_roi_repeat|	perform whole genome mapping + slice on ROI + longshot phasing + repeat calling STRique|--strique_config \<strique_config\> --roi \<roi\>|
|wgs_splitcas9_repeat|	perform whole genome mapping + split based on Cas9 sites + longshot phasing + repeat calling STRique|--splitfile \<splitfile\> --strique_config \<strique_config\>|
|SMA_splitcas9|perform whole genome mapping + split based on Cas9 sites +  SMA Variant calling + haplotype phasing|--roi \<roi\> --splitfile \<splitfileSMA\>|
|SMA_adaptive|perform whole genome mapping (if required) +  SMA Variant calling + haplotype phasing|--ploidy \<SMN2 copy number\>|

in which:

| Method | Description | Format | 
| :--- | :--- | :--- |
|\<strique_config\>|STRique config file|see STRIque documentation for correct format|
|\<roi\>|Region of Interest|\<chromosome\>:\<from\>-\<to\>|
|\<splitfile\>|Tab delimited file containing Cas9 split sites for each sample (may contain more lines per sample)|\<SampleID\> \<chromosome\> \<position1\> \<postion2\>|
|\<splitfileSMA\>|Tab delimited file containing Cas9 split sites for each sample (may contain more lines per sample) and SMN2 copynumber (ploidy)|\<SampleID\> \<chromosome> \<position1\> \<postion2\> \<ploidy\>|
|\<SMN2 copy number\>|SMN2 copynumber (ploidy)|\<int\>|

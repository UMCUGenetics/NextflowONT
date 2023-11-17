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
bam_single = start from single bam file (should be only bam in the folder) without any Guppy/ONT information.
\

| Method | Description | Optional parameters needed|
| --- | :--- | :--- |
|wgs|perform whole genome mapping + longshot phasing|none|
|SMA_adaptive|perform whole genome mapping (if required) +  SMA Variant calling + haplotype phasing|--ploidy \<SMN2 copy number\>|

in which:

| Method | Description | Format | 
| :--- | :--- | :--- |
|\<SMN2 copy number\>|SMN2 copynumber (ploidy)|\<int\>|

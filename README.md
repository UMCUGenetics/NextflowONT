# NextflowONT
Nextflow Oxford Nanopore Technologies workflow
Note this workflow is tested on R9.4.1 sequencing data with guppy_6.1.2 basecalling only.

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
nextflow run ONT.nf -c ONT.config --input_path <input_path> --outdir <output_dir_path> --start <bam|rebase> --method <method> --email <email> [-profile slurm|mac]
```
###<input_path>
Input path is full path to either FAST5 or BAM file, depending on method (see below)
Rebase assumes raw data folder with fast5_pass and fast5_fail subfolder.
bam/bam_remap assumes Guppy output folder with subfolder pass in which BAM files are located
bam_single/bam_single_remap assumes a single BAM file per input folder
###<methods>
rebase = include re-basecalling\
bam = start from Guppy folder including bam files.\
bam_remap = start from Guppy folder including bam files, but perform remapping to genome in config.
bam_single = start from single bam file (should be only bam in the folder) without any Guppy/ONT information.
bam_single_remap = start from single bam file (should be only bam in the folder) without any Guppy/ONT information and perform remapping.

| Method | Description | Optional parameters needed|
| --- | :--- | :--- |
|wgs|perform whole genome mapping + longshot phasing|none|
|SMA_adaptive|perform whole genome mapping (if required) +  SMA Variant calling + haplotype phasing|--ploidy \<SMN2 copy number\>|

in which:

| Method | Description | Format | 
| :--- | :--- | :--- |
|\<SMN2 copy number\>|SMN2 copynumber (ploidy)|\<int\>|

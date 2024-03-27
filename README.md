# NextflowONT
Nextflow Oxford Nanopore Technologies workflow\
Note: this workflow hase been tested on R9.4.1 sequencing data with guppy_6.1.2 basecalling only.\
Note: current workflow setup has been tested on Rocky Linux 8 combined with Slurm Workload Manager.

#### Configure paths before use
ONT_wgs_all.config:
<pre>
  genome_fasta          full path to reference genome fasta (.fasta/.fa/.fna)
  genome_mapping_index  full path to reference genome minimap2 index (.mmi)
  calling_target_bed    full path to "position specific" variant calling (.bed)
  calling_target_region region interest for variant calling based on reference genome (chr:start-stop)
  homopolymer_bed       full path to homopolymer region of reference genome (.bed)
  guppy_basecaller_path full path to guppy_basecaller executable
  guppy_path            full path to guppy folder
  runOptions            run options for singularity
  cacheDir              singularity image cache folder
</pre>

#### Running ONT workflow
```bash
nextflow run ONT.nf -c ONT.config --input_path <input_path> --outdir <output_dir_path> --start <bam|rebase> --method <method> --email <email> [-profile slurm]
```

##### --input_path
<pre>
Input path is full path to either FAST5 or BAM file, depending on method (see below).
  * rebase              assumes raw data folder with fast5_pass and fast5_fail subfolder.
  * bam(_remap)         assumes Guppy output folder with subfolder pass in which BAM files are located.
  * bam_single(_remap)  assumes a single BAM file per input folder.
</pre>
##### --method
<pre>
  * rebase              start from FAST5 raw data ouput folder and include re-basecalling.
  * bam                 start from Guppy folder including bam files.
  * bam_remap           start from Guppy folder including bam files, but perform remapping to genome in config.
  * bam_single          start from single bam file (should be only bam in the folder) without any Guppy/ONT information.
  * bam_single_remap    start from single bam file (should be only bam in the folder) without any Guppy/ONT information and perform remapping.
</pre>

| Method | Description | Optional parameters needed|
| --- | :--- | :--- |
|wgs|perform whole genome mapping + longshot phasing|none|
|SMA_adaptive|perform whole genome mapping (if required) +  SMA Variant calling + haplotype phasing|--ploidy \<SMN2 copy number\>|

in which:

| Method | Description | Format |
| :--- | :--- | :--- |
|\<SMN2 copy number\>|SMN2 copynumber (ploidy)|\<int\>|

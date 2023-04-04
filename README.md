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
rebase = include re-basecalling
bam = start from Guppy folder indcluding bam files.
<method> = 
	wgs			perform whole genome mapping + longshot phasing.
	wgs_repeat		perform whole genome mapping + longshot phasing + repeat calling STRique. Needs optional parameter --strique_config <strique_config>
	wgs_roi			perform whole genome mapping + slice on ROI + longshot phasing.  Needs optional parameter --roi <roi_file>
	wgs_roi_repeat		perform whole genome mapping + slice on ROI + longshot phasing + repeat calling STRique. Needs optional parameters --strique_config <strique_config> --roi <roi_file>
	wgs_splitcas9_repeat 	perform whole genome mapping + split based on Cas9 sites + longshot phasing + repeat calling STRique. Needs optional parameters --splitfile <splitfile_repeat> --strique_config <strique_config>
	targeted		perform whole genome mapping + targeted remapping (SMN2 default in .config). Needs optional paramete --roi <roi_file>
	targeted_splitcas9	perform whole genome mapping + targeted remapping (SMN2 default in .config) + split BAM based on Cas9 sites. Needs optional parameters --roi <roi file> --splitfile <splitfile>
        targeted_SMA_splitcas9	perform whole genome mapping + split based on Cas9 sites +  SMA Variant calling + haplotype phasing. Needs optional parameters --roi <roi> --splitfile <splitfileSMA>
	targeted_SMA_adaptive	perform whole genome mapping + targeted remapping (SMN2 default in .config) +  SMA Variant calling + haplotype phasing. Needs optional parameters --roi <roi> --ploidy <SMN2 copy number>

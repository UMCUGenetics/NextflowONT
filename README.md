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
nextflow run ONT.nf -c ONT.config --fastq_path <fastq_dir_path> --outdir <output_dir_path> --email <email> --roi <BED file of ROI> [-profile slurm|mac]
```


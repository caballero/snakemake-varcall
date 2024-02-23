# Simple Snakemake pipeline to do Variant Calling

## Process

![workflow](img/workflow.svg)

## Preparations

1. Genome needs to be indexed for BWA, fasta and index is defined in config file
2. Fastq files are located in a directory defined in config, it expects a `SAMPLEID_1.fastq.gz` and `SAMPLEID_2.fastq.gz` files
3. All processes and results are defined by `-d` parameter in snakemake
4. All logs are stored in `log/` directory
5. Dependencies can be added in a conda environment as defined in `environment.yaml` for running in local mode, or defined individually in `envs` to use with `--use-conda`or as modules to be used with `--use-envmodules`
6. The pipeline can be run in a cluster (`--slurm`) if configured properly

## Input configuration

1. Defined in config.yaml, it requires

```
samples:
  - sampleID1
  - sampleID2
fastq_dir: /path/to/fastq/dir
genome: /path/to/genome.fa
bcf_mpileup_param: "params for bcftools mpileup"
bcf_call_param: "params for bcftools mpileup"
bcf_norm_param: "params for bcftools norm"
bcf_filter_param: "params for bcftools filters"
```

## Execution

`snakemake <params>`

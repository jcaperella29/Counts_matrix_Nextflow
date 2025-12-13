Counts_matrix_Nextflow

A minimal, portable Nextflow DSL2 RNA-seq pipeline that:

Accepts paired-end FASTQs

Performs QC → alignment → counting → matrix merge

Works with Conda, Docker, or Singularity

Uses a toy reference FASTA by default

Automatically builds BWA indexes if missing

Allows users to provide their own reference genome

This is intentionally lightweight and ideal for:

Teaching / demos

Testing infrastructure

Prototyping larger RNA-seq workflows

Pipeline Overview

Steps:

(Optional) Generate toy reference FASTA (MAKE_REF)

Build BWA index if missing (INDEX_REF)

Paired-end QC with fastp

Paired-end alignment with BWA MEM

Sorting & indexing with samtools

Counting with featureCounts (SAF mode)

Merge per-sample counts into a matrix


data/
├── sample_1/
│   ├── sample_1_R1.fastq.gz
│   └── sample_1_R2.fastq.gz
└── sample_2/
    ├── sample_2_R1.fastq.gz
    └── sample_2_R2.fastq.gz

Supported FASTQ name patterns:
Supported FASTQ name patterns:

*_R1*.fastq.gz / *_R2*.fastq.gz

*1.fastq.gz / *2.fastq.gz

Outputs
results/
├── ref/
│   ├── ref.fa
│   └── ref.fa.*        # BWA index files
├── qc/
│   ├── sample_1.fastp.html
│   └── sample_1.fastp.json
├── bam/
│   ├── sample_1.bam
│   └── sample_1.bam.bai
├── counts_per_sample/
│   ├── sample_1.counts.tsv
│   └── sample_2.counts.tsv
└── counts_matrix.tsv


Running the Pipeline(all the code is in bash)
1️⃣ Conda (recommended for local dev)
conda activate nf_genomics
nextflow run main.nf -profile conda


2️⃣ Docker

Make sure Docker is running and the image exists locally:

docker build -t rnaseq-pipeline .
nextflow run main.nf -profile docker

3️⃣ Singularity / Apptainer (HPC-friendly)

sudo singularity build containers/rnaseq-pipeline.sif docker-archive://containers/rnaseq-pipeline.tar
nextflow run main.nf -profile singularity

Using Your Own Reference Genome

By default, the pipeline generates a toy FASTA (ref.fa).

To use a real reference:

nextflow run main.nf \
  -profile conda \
  --ref path/to/genome.fa

BWA index files will be created automatically if missing

Existing indexes will be reused


Parameters

| Parameter     | Default   | Description              |
| ------------- | --------- | ------------------------ |
| `--input_dir` | `data`    | Sample directories       |
| `--outdir`    | `results` | Output directory         |
| `--threads`   | `4`       | Threads per task         |
| `--ref`       | *(none)*  | Optional reference FASTA |



Profiles Summary


| Profile       | Use case                |
| ------------- | ----------------------- |
| `conda`       | Local development       |
| `docker`      | Reproducible containers |
| `singularity` | HPC clusters            |

Configured in nextflow.config.

equirements

Nextflow ≥ 23

One of:

Conda

Docker

Singularity / Apptainer

Notes & Design Philosophy

This pipeline is intentionally simple

SAF annotation is generated automatically for toy refs

FeatureCounts is run in paired-end mode

Ideal as a foundation for:

GTF-based workflows

Salmon / STAR integration

Multi-omics pipelines


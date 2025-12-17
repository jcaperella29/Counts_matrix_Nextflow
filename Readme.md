Counts_matrix_Nextflow

A portable Nextflow DSL2 RNA-seq pipeline that performs:

QC

Genome alignment

Gene-level counting

Transcript-level quantification

Coverage track generation

Matrix merging and summary reporting

Designed to be minimal, readable, and robust, while still demonstrating real RNA-seq best practices.

Features

This pipeline:

Accepts paired-end FASTQs via a samplesheet

Performs QC → alignment → quantification → matrix merge

Generates gene counts, TPM matrices, and BigWig tracks

Automatically builds STAR and Salmon indexes

Skips coverage generation gracefully for samples with zero mapped reads

Works with Conda, Docker, or Singularity/Apptainer

Uses standard, well-documented tools

Ideal for:

Teaching / demos

Testing infrastructure

Prototyping larger RNA-seq workflows

Portfolio or template pipelines

Pipeline Overview
Steps

FASTQ QC & trimming

fastp

Genome index build (once per run)

STAR --runMode genomeGenerate

Genome alignment

STAR → coordinate-sorted BAM

BAM indexing + flagstat

Coverage tracks

deepTools bamCoverage → BigWig

Automatically skipped if no mapped reads

Gene-level counting

featureCounts (GTF-based)

Transcript-level quantification

Salmon quant

Matrix merging

Gene count matrix (counts_matrix.tsv)

Transcript TPM matrix (salmon_tpm_matrix.tsv)

Summary reporting

MultiQC

Input Format
Samplesheet (CSV)
sample,read1,read2
S1,data/S1_R1.fastq.gz,data/S1_R2.fastq.gz
S2,data/S2_R1.fastq.gz,data/S2_R2.fastq.gz


Required columns:

sample

read1

read2

Required Inputs

Reference genome FASTA (--ref)

Gene annotation GTF (--gtf)

Transcript FASTA for Salmon (--transcripts)

Samplesheet CSV (--samplesheet)

Output Structure
results/
├── qc/
│   ├── S1.fastp.html
│   └── S1.fastp.json
├── ref/
│   ├── star/
│   │   └── STAR_INDEX/
│   └── salmon/
│       └── SALMON_INDEX/
├── bam/
│   ├── S1.bam
│   ├── S1.bam.bai
│   └── S1.flagstat.txt
├── bigwig/
│   └── S1.bw
├── counts_per_sample/
│   ├── S1.counts.tsv
│   └── S2.counts.tsv
├── counts_matrix.tsv
├── salmon_tpm_matrix.tsv
└── multiqc_report.html

Deployment options (all code is bash)
 to run  with Docker
 docker build -t rnaseq-pipeline .
nextflow run main.nf -profile docker \
  --samplesheet samples.csv \
  --ref genome.fa \
  --gtf genes.gtf \
  --transcripts transcripts.fa

to run with Apptainer

singularity build containers/rnaseq-pipeline.sif docker://rnaseq-pipeline
nextflow run main.nf -profile singularity \
  --samplesheet samples.csv \
  --ref genome.fa \
  --gtf genes.gtf \
  --transcripts transcripts.fa
Parameters
| Parameter       | Default      | Description                   |
| --------------- | ------------ | ----------------------------- |
| `--samplesheet` | *(required)* | CSV mapping samples to FASTQs |
| `--ref`         | *(required)* | Reference genome FASTA        |
| `--gtf`         | *(required)* | Gene annotation GTF           |
| `--transcripts` | *(required)* | Transcript FASTA for Salmon   |
| `--outdir`      | `results`    | Output directory              |
| `--threads`     | `4`          | Threads per task              |
| `--bw_binsize`  | `10`         | BigWig bin size               |
| `--bw_norm`     | `CPM`        | BigWig normalization          |

Requirements

Nextflow ≥ 23

Docker 
Singularity / Apptainer
fastp
STAR
samtools
subread (featureCounts)
salmon
deeptools
multiqc
otes & Design Philosophy

Intentionally simple and readable

Avoids over-engineering

Uses explicit channels instead of heavy abstraction

Guards against common RNA-seq failure modes (e.g. zero-mapped samples)

Designed to be:

Extended with DESeq2 / edgeR

Modularized into modules/

Integrated into larger multi-omics workflows




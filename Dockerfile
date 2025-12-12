# Dockerfile for mini RNA-seq Nextflow pipeline
# Provides: python, fastp, bwa, samtools, subread (featureCounts), gzip

FROM condaforge/mambaforge:latest

# Create a dedicated env from the same spec as envs/rnaseq.yml
# (Assumes envs/rnaseq.yml is copied into the image by the build context.)

WORKDIR /opt/pipeline

COPY envs/rnaseq.yml /opt/pipeline/envs/rnaseq.yml

RUN mamba env create -f /opt/pipeline/envs/rnaseq.yml && \
    mamba clean -afy

# Put the rnaseq env on PATH so Nextflow processes see the tools directly
ENV PATH="/opt/conda/envs/rnaseq/bin:${PATH}"

# Optional: set a default workdir inside the container
WORKDIR /workspace

# Default command (Nextflow overrides this when running processes)
CMD ["bash"]

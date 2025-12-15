# Dockerfile for portfolio-ready RNA-seq Nextflow pipeline
# Tools: python, fastp, STAR, samtools, subread(featureCounts), salmon, multiqc, deeptools

FROM condaforge/mambaforge:latest

WORKDIR /opt/pipeline

COPY envs/rnaseq.yml /opt/pipeline/envs/rnaseq.yml

# Build env from your YAML, then install the additional tools needed for
# STAR+Salmon+BigWig tracks+MultiQC.
RUN mamba env create -f /opt/pipeline/envs/rnaseq.yml && \
    mamba install -n rnaseq -c bioconda -c conda-forge \
      star \
      salmon \
      multiqc \
      deeptools \
      && \
    mamba clean -afy

# Put the rnaseq env on PATH so Nextflow processes see the tools directly
ENV PATH="/opt/conda/envs/rnaseq/bin:${PATH}"

WORKDIR /workspace

CMD ["bash"]

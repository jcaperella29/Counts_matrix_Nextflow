#!/usr/bin/env bash
#SBATCH --job-name=nf_rnaseq
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=general

# If your cluster requires a specific account:
##SBATCH --account=YOUR_ACCOUNT

set -euo pipefail

# --- User-editable paths ---
WORKDIR=${WORKDIR:-$PWD}
INPUT_DIR=${INPUT_DIR:-data}
REF=${REF:-ref.fa}
OUTDIR=${OUTDIR:-results}
THREADS=${THREADS:-8}
PROFILE=${PROFILE:-singularity}

# Where Nextflow should put its cache + work directory (cluster-friendly)
NXF_HOME=${NXF_HOME:-$WORKDIR/.nextflow}
NXF_WORK=${NXF_WORK:-$WORKDIR/work}

# Optional: where Singularity caches images
SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-$WORKDIR/.singularity_cache}

mkdir -p "$WORKDIR/logs" "$NXF_HOME" "$NXF_WORK" "$SINGULARITY_CACHEDIR"

cd "$WORKDIR"

# --- Load modules (edit for your cluster) ---
# module purge
# module load nextflow
# module load singularity

export NXF_HOME
export NXF_WORK
export NXF_OPTS='-Xms1g -Xmx4g'
export SINGULARITY_CACHEDIR

# Recommended for SLURM + Nextflow
export NXF_ANSI_LOG=false

# Run pipeline
nextflow run main.nf \
  -profile "$PROFILE" \
  --input_dir "$INPUT_DIR" \
  --ref "$REF" \
  --outdir "$OUTDIR" \
  --threads "$THREADS" \
  -with-report "$OUTDIR/nf_report.html" \
  -with-trace  "$OUTDIR/nf_trace.txt" \
  -with-timeline "$OUTDIR/nf_timeline.html"

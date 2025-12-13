nextflow.enable.dsl = 2

params.input_dir = params.input_dir ?: 'data'
params.outdir    = params.outdir    ?: 'results'
params.threads   = params.threads   ?: 4
params.ref       = params.ref ?: null



process INDEX_REF {

    tag "bwa_index"
    publishDir "${params.outdir}/ref", mode: 'copy'

    input:
    path ref

    output:
    tuple path(ref),
          path("${ref}.bwt"),
          path("${ref}.sa"),
          path("${ref}.ann"),
          path("${ref}.amb"),
          path("${ref}.pac")

    script:
    """
    echo "Checking BWA index for ${ref}"

    if [ ! -f "${ref}.bwt" ]; then
        echo "BWA index not found — building index"
        bwa index ${ref}
    else
        echo "BWA index already exists — skipping"
    fi
    """
}
/* -------------------------------------------------------------------------- */
/* MAKE_REF – generate a minimal reference FASTA for toy/testing runs          */
/* -------------------------------------------------------------------------- */
process MAKE_REF {

    tag "ref"
    publishDir "${params.outdir}/ref", mode: 'copy'

    output:
    path "ref.fa"

    script:
    """
    cat > ref.fa << 'EOF'
>chr1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
EOF
    """
}


/* -------------------------------------------------------------------------- */
/* FASTP_QC – per-sample fastp QC                                             */
/* -------------------------------------------------------------------------- */
process FASTP_QC {

    tag "${sample_id}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(sample_dir)

    output:
    path "${sample_id}.fastp.html"
    path "${sample_id}.fastp.json"

    script:
    """
    echo "Running fastp QC for sample: ${sample_id}"

    R1=\$(ls ${sample_dir}/*_R1*.fastq.gz ${sample_dir}/*1.fastq.gz 2>/dev/null | head -n1)
    R2=\$(ls ${sample_dir}/*_R2*.fastq.gz ${sample_dir}/*2.fastq.gz 2>/dev/null | head -n1)

    if [ -z "\$R1" ] || [ -z "\$R2" ]; then
        echo "ERROR: Could not find R1/R2 FASTQs for ${sample_id}" >&2
        exit 1
    fi

    fastp \\
      -i \$R1 -I \$R2 \\
      -h ${sample_id}.fastp.html \\
      -j ${sample_id}.fastp.json \\
      -w ${params.threads}
    """
}


/* -------------------------------------------------------------------------- */
/* ALIGN_READS – BWA paired-end alignment + sorting + BAM index               */
/* -------------------------------------------------------------------------- */
process ALIGN_READS {

    tag "${sample_id}"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(sample_dir), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    echo "Aligning sample: ${sample_id}"
    echo "Using reference: ${ref}"

    R1=\$(ls ${sample_dir}/*_R1*.fastq.gz ${sample_dir}/*1.fastq.gz 2>/dev/null | head -n1)
    R2=\$(ls ${sample_dir}/*_R2*.fastq.gz ${sample_dir}/*2.fastq.gz 2>/dev/null | head -n1)

    if [ -z "\$R1" ] || [ -z "\$R2" ]; then
        echo "ERROR: Missing FASTQs for ${sample_id}" >&2
        exit 1
    fi

    echo "R1: \$R1"
    echo "R2: \$R2"

    # Build BWA index in task work dir (toy pipeline)
    bwa index ${ref}

    # Align and sort
    bwa mem -t ${params.threads} ${ref} \$R1 \$R2 \\
      | samtools sort -@ ${params.threads} -o ${sample_id}.bam

    samtools index ${sample_id}.bam
    """
}


/* -------------------------------------------------------------------------- */
/* FEATURECOUNTS_PER_SAMPLE – paired-end featureCounts with auto-SAF          */
/* -------------------------------------------------------------------------- */
process FEATURECOUNTS_PER_SAMPLE {

    tag "${sample_id}"
    publishDir "${params.outdir}/counts_per_sample", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.counts.tsv")

    script:
    """
    echo "Running featureCounts for sample: ${sample_id}"
    echo "BAM: ${bam}"

    # Minimal SAF annotation for toy reference
    cat > annotation.saf << 'EOF'
GeneID\tChr\tStart\tEnd\tStrand
gene1\tchr1\t1\t1000000\t+
EOF

    featureCounts \\
      -T ${params.threads} \\
      -p \\
      -B \\
      -F SAF \\
      -a annotation.saf \\
      -o ${sample_id}.raw_counts.txt \\
      ${bam}

    # Extract (feature, count) into sample_X.counts.tsv
    awk 'BEGIN{OFS="\\t"} \\
         /^#/ {next} \\
         NR==2 {next} \\
         NR>2 {print \$1, \$7}' ${sample_id}.raw_counts.txt \\
      > ${sample_id}.counts.tsv
    """
}


/* -------------------------------------------------------------------------- */
/* MERGE_COUNTS_MATRIX – merge sample-level TSVs into one matrix              */
/* -------------------------------------------------------------------------- */
process MERGE_COUNTS_MATRIX {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(count_files)

    output:
    path "counts_matrix.tsv"

    script:
    """
    python << 'EOF'
import glob
import os

files = sorted(glob.glob("*.counts.tsv"))
if not files:
    raise SystemExit("No count files found.")

all_counts = {}
features = set()

for fname in files:
    sid = os.path.basename(fname).split(".")[0]
    all_counts[sid] = {}
    with open(fname) as f:
        for line in f:
            feat, count = line.strip().split("\\t")
            all_counts[sid][feat] = count
            features.add(feat)

features = sorted(features)

with open("counts_matrix.tsv", "w") as out:
    out.write("feature\\t" + "\\t".join(all_counts.keys()) + "\\n")
    for feat in features:
        row = [feat] + [all_counts[sid].get(feat, "0") for sid in all_counts]
        out.write("\\t".join(row) + "\\n")
EOF
    """
}


/* -------------------------------------------------------------------------- */
workflow {

    file(params.outdir).mkdirs()

    // Decide reference source (toy by default, user can override with --ref)
    if (params.ref) {
        log.info "Using user-provided reference: ${params.ref}"
        raw_ref_ch = Channel.fromPath(params.ref, checkIfExists: true)
    } else {
        log.info "No reference provided -- using toy reference"
        raw_ref_ch = MAKE_REF()
    }

    // Build BWA index if missing (runs once per resume-able execution)
    indexed_ref_ch = raw_ref_ch | INDEX_REF

    // We only need to pass the ref fasta path into ALIGN_READS
    ref_only_ch = indexed_ref_ch.map { ref, bwt, sa, ann, amb, pac -> ref }

    // sample dirs: (sample_id, sample_dir)
    Channel
        .fromPath("${params.input_dir}/*", type: 'dir')
        .map { dir -> tuple(dir.baseName, dir) }
        .set { samples_ch }

    // QC (branch)
    samples_ch | FASTP_QC

    // Align
    aligned_ch = samples_ch
        .combine(ref_only_ch)
        | ALIGN_READS

    // Count per sample
    counts_ch = aligned_ch | FEATURECOUNTS_PER_SAMPLE

    // Merge counts matrix
    counts_ch
        .map { sid, f -> f }
        .collect()
        .set { counts_files_ch }

    counts_files_ch | MERGE_COUNTS_MATRIX
}


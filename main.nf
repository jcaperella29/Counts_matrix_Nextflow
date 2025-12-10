nextflow.enable.dsl = 2

params.input_dir = params.input_dir ?: 'data'
params.outdir    = params.outdir    ?: 'results'
params.ref       = params.ref       ?: 'ref.fa'   // reference genome
params.threads   = params.threads   ?: 4


/*
 * Process: ALIGN_READS
 * Input:  (sample_id, sample_dir, ref)
 * Output: (sample_id, sample_id.bam)
 */
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

    # Try to find read pairs; adapt these patterns as needed.
    R1=\$(ls ${sample_dir}/*_R1*.fastq.gz ${sample_dir}/*1.fastq.gz 2>/dev/null | head -n1)
    R2=\$(ls ${sample_dir}/*_R2*.fastq.gz ${sample_dir}/*2.fastq.gz 2>/dev/null | head -n1)

    if [ -z "\$R1" ] || [ -z "\$R2" ]; then
        echo "Could not find R1/R2 FASTQs for ${sample_id} in ${sample_dir}" >&2
        exit 1
    fi

    echo "R1: \$R1"
    echo "R2: \$R2"

    # Build BWA index for this reference in the work dir
    bwa index ${ref}

    # Align and sort
    bwa mem -t ${params.threads} ${ref} \$R1 \$R2 \\
      | samtools sort -@ ${params.threads} -o ${sample_id}.bam

    samtools index ${sample_id}.bam
    """
}


/*
 * Process: COUNT_FROM_BAM
 * Input:  (sample_id, bam)
 * Output: (sample_id, sample_id.counts.tsv)
 */
process COUNT_FROM_BAM {

    tag "${sample_id}"
    publishDir "${params.outdir}/counts_per_sample", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.counts.tsv")

    script:
    """
    echo "Counting reads for sample: ${sample_id}"

    samtools idxstats ${bam} \\
      | awk 'BEGIN{OFS="\\t"} {print \$1, \$3}' \\
      > ${sample_id}.counts.tsv
    """
}


/*
 * Process: MERGE_COUNTS_MATRIX
 * Input:  list of count_files (Nextflow stages them all here)
 * Output: counts_matrix.tsv
 *
 * Uses ONLY built-in Python, no pandas.
 */
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

# Find all staged count files
files = sorted(glob.glob("*.counts.tsv"))
if not files:
    raise SystemExit("No *.counts.tsv files found")

# Load all counts into dictionaries
all_counts = {}
features = set()

for fname in files:
    sid = os.path.basename(fname).split(".")[0]
    all_counts[sid] = {}
    with open(fname) as fh:
        for line in fh:
            feat, count = line.strip().split("\\t")
            all_counts[sid][feat] = count
            features.add(feat)

features = sorted(list(features))

# Write merged matrix
with open("counts_matrix.tsv", "w") as out:
    header = ["feature"] + list(all_counts.keys())
    out.write("\\t".join(header) + "\\n")

    for feat in features:
        row = [feat]
        for sid in all_counts.keys():
            row.append(all_counts[sid].get(feat, "0"))
        out.write("\\t".join(row) + "\\n")
EOF
    """
}


/*
 * Main workflow
 */
workflow {

    // Ensure output directory exists
    file(params.outdir).mkdirs()

    /*
     * Channel of samples:
     *   (sample_id, sample_dir)
     */
    Channel
        .fromPath("${params.input_dir}/*", type: 'dir')
        .map { dir ->
            tuple(dir.baseName, dir)
        }
        .set { samples_ch }

    /*
     * Single-value channel with the reference FASTA
     */
    Channel
        .value( file(params.ref) )
        .set { ref_ch }

    /*
     * Combine samples with reference:
     *   (sample_id, sample_dir, ref)
     */
    samples_with_ref_ch = samples_ch.combine(ref_ch)

    /*
     * 1) Align → (sample_id, bam)
     */
    aligned_ch = samples_with_ref_ch | ALIGN_READS

    /*
     * 2) Count → (sample_id, sample_id.counts.tsv)
     */
    counts_ch  = aligned_ch | COUNT_FROM_BAM

    /*
     * 3) Collect all count files into a single list channel
     *    and feed to MERGE_COUNTS_MATRIX
     */
    counts_ch
        .map { sid, f -> f }    // keep only the file path
        .collect()              // list of paths, single emission
        .set { counts_files_ch }

    counts_files_ch | MERGE_COUNTS_MATRIX
}

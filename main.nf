nextflow.enable.dsl = 2

/*
  Portfolio-ready mini RNA-seq pipeline
  - QC: fastp
  - Genome alignment (for tracks + gene counts): STAR -> sorted BAM
  - Tracks: BigWig (deeptools bamCoverage)
  - Gene counts: featureCounts (GTF)
  - Transcript quant: Salmon
  - Summary: MultiQC

  Required tools in container:
    fastp, star, samtools, subread(featureCounts), salmon, multiqc, deeptools

  Inputs (recommended): samplesheet.csv with columns: sample,read1,read2
*/

params.samplesheet = params.samplesheet ?: null
params.input_dir   = params.input_dir   ?: 'data'     // fallback if no samplesheet
params.outdir      = params.outdir      ?: 'results'
params.threads     = params.threads     ?: 4

// Reference inputs
params.ref         = params.ref         ?: null        // genome FASTA
params.gtf         = params.gtf         ?: null        // annotation GTF for featureCounts
params.transcripts = params.transcripts ?: null        // transcriptome FASTA for Salmon

// Track settings
params.bw_binsize  = params.bw_binsize  ?: 10
params.bw_norm     = params.bw_norm     ?: 'CPM'       // CPM is a reasonable default

/* -------------------------------------------------------------------------- */
/* Helpers                                                                    */
/* -------------------------------------------------------------------------- */

def require_param(name, value) {
    if (!value) {
        throw new IllegalArgumentException("Missing required parameter --${name}")
    }
}

/* -------------------------------------------------------------------------- */
/* MAKE_REF – generate a minimal reference FASTA for toy/testing runs          */
/* (kept for completeness; real runs should provide --ref/--gtf/--transcripts) */
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
    publishDir "${params.outdir}/qc/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trim.fastq.gz"),
          path("${sample_id}_R2.trim.fastq.gz"),
          path("${sample_id}.fastp.html"),
          path("${sample_id}.fastp.json")

    script:
    """
    echo "Running fastp for: ${sample_id}"

    fastp \\
      -i ${r1} -I ${r2} \\
      -o ${sample_id}_R1.trim.fastq.gz \\
      -O ${sample_id}_R2.trim.fastq.gz \\
      -h ${sample_id}.fastp.html \\
      -j ${sample_id}.fastp.json \\
      -w ${params.threads}
    """
}

/* -------------------------------------------------------------------------- */
/* STAR_INDEX – build STAR genome index                                       */
/* -------------------------------------------------------------------------- */
process STAR_INDEX {

    tag "star_index"
    publishDir "${params.outdir}/ref/star", mode: 'copy'

    input:
    path ref

    output:
    path "STAR_INDEX"

    script:
    """
    mkdir -p STAR_INDEX
    STAR \\
      --runThreadN ${params.threads} \\
      --runMode genomeGenerate \\
      --genomeDir STAR_INDEX \\
      --genomeFastaFiles ${ref}
    """
}

/* -------------------------------------------------------------------------- */
/* STAR_ALIGN – align paired-end reads, output coordinate-sorted BAM          */
/* -------------------------------------------------------------------------- */
process STAR_ALIGN {

    tag "${sample_id}"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2), path(star_index_dir)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    echo "STAR aligning: ${sample_id}"

    STAR \\
      --runThreadN ${params.threads} \\
      --genomeDir ${star_index_dir} \\
      --readFilesIn ${r1} ${r2} \\
      --readFilesCommand zcat \\
      --outSAMtype BAM SortedByCoordinate \\
      --outFileNamePrefix ${sample_id}.

    # STAR writes: <prefix>Aligned.sortedByCoord.out.bam
    mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.bam
    samtools index ${sample_id}.bam

    # Basic mapping stats for MultiQC
    samtools flagstat ${sample_id}.bam > ${sample_id}.flagstat.txt
    """
}

/* -------------------------------------------------------------------------- */
/* BAM_TO_BIGWIG – genome browser tracks via deeptools bamCoverage            */
/* -------------------------------------------------------------------------- */
process BAM_TO_BIGWIG {

    tag "${sample_id}"
    publishDir "${params.outdir}/tracks", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bw")

    script:
    """
    echo "Creating BigWig track for: ${sample_id}"

    bamCoverage \\
      -b ${bam} \\
      -o ${sample_id}.bw \\
      --binSize ${params.bw_binsize} \\
      --normalizeUsing ${params.bw_norm} \\
      -p ${params.threads}
    """
}

/* -------------------------------------------------------------------------- */
/* FEATURECOUNTS_PER_SAMPLE – gene-level counts from STAR BAM                 */
/* -------------------------------------------------------------------------- */
process FEATURECOUNTS_PER_SAMPLE {

    tag "${sample_id}"
    publishDir "${params.outdir}/counts_per_sample", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(gtf)

    output:
    tuple val(sample_id), path("${sample_id}.counts.tsv"), path("${sample_id}.featureCounts.summary")

    script:
    """
    echo "Running featureCounts for: ${sample_id}"

    featureCounts \\
      -T ${params.threads} \\
      -p \\
      -B \\
      -a ${gtf} \\
      -o ${sample_id}.raw_counts.txt \\
      ${bam}

    # Keep the summary too (MultiQC can parse it)
    cp ${sample_id}.raw_counts.txt.summary ${sample_id}.featureCounts.summary

    # Extract (feature, count)
    awk 'BEGIN{OFS="\t"} \\
         /^#/ {next} \\
         NR==2 {next} \\
         NR>2 {print $1, $7}' ${sample_id}.raw_counts.txt \\
      > ${sample_id}.counts.tsv
    """
}

/* -------------------------------------------------------------------------- */
/* MERGE_COUNTS_MATRIX – merge sample gene count TSVs                         */
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
            feat, count = line.strip().split("\t")
            all_counts[sid][feat] = count
            features.add(feat)

features = sorted(features)

with open("counts_matrix.tsv", "w") as out:
    out.write("feature\t" + "\t".join(all_counts.keys()) + "\n")
    for feat in features:
        row = [feat] + [all_counts[sid].get(feat, "0") for sid in all_counts]
        out.write("\t".join(row) + "\n")
EOF
    """
}

/* -------------------------------------------------------------------------- */
/* SALMON_INDEX – build Salmon transcriptome index                            */
/* -------------------------------------------------------------------------- */
process SALMON_INDEX {

    tag "salmon_index"
    publishDir "${params.outdir}/ref/salmon", mode: 'copy'

    input:
    path transcripts

    output:
    path "SALMON_INDEX"

    script:
    """
    mkdir -p SALMON_INDEX
    salmon index \\
      -t ${transcripts} \\
      -i SALMON_INDEX \\
      -p ${params.threads}
    """
}

/* -------------------------------------------------------------------------- */
/* SALMON_QUANT – per-sample Salmon quant                                     */
/* -------------------------------------------------------------------------- */
process SALMON_QUANT {

    tag "${sample_id}"
    publishDir "${params.outdir}/salmon/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2), path(salmon_index)

    output:
    tuple val(sample_id), path("quant.sf")

    script:
    """
    echo "Salmon quant for: ${sample_id}"

    salmon quant \\
      -i ${salmon_index} \\
      -l A \\
      -1 ${r1} -2 ${r2} \\
      -p ${params.threads} \\
      --validateMappings \\
      -o .

    # keep a stable filename for downstream merge
    cp quant.sf quant.sf
    """
}

/* -------------------------------------------------------------------------- */
/* MERGE_SALMON_TPM – merge Salmon quant.sf into TPM matrix                   */
/* -------------------------------------------------------------------------- */
process MERGE_SALMON_TPM {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(qsf)

    output:
    path "salmon_tpm_matrix.tsv"

    script:
    """
    python << 'EOF'
import os
import glob

# Expect files placed in workdir with names like <sample>.quant.sf via staging
qfiles = sorted(glob.glob("*.quant.sf"))
if not qfiles:
    raise SystemExit("No *.quant.sf files found for merging")

samples = []
rows = {}  # tx -> {sample: tpm}

for q in qfiles:
    sample = os.path.basename(q).replace(".quant.sf", "")
    samples.append(sample)
    with open(q) as f:
        header = f.readline().strip().split("\t")
        # Name Length EffectiveLength TPM NumReads
        idx_name = header.index("Name")
        idx_tpm  = header.index("TPM")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            tx = parts[idx_name]
            tpm = parts[idx_tpm]
            rows.setdefault(tx, {})[sample] = tpm

samples = sorted(samples)

with open("salmon_tpm_matrix.tsv", "w") as out:
    out.write("transcript\t" + "\t".join(samples) + "\n")
    for tx in sorted(rows.keys()):
        out.write(tx)
        for s in samples:
            out.write("\t" + rows[tx].get(s, "0"))
        out.write("\n")
EOF
    """
}

/* -------------------------------------------------------------------------- */
/* MULTIQC – aggregate QC/alignment/counting reports                          */
/* -------------------------------------------------------------------------- */
process MULTIQC {

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path(qc_dir)

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc ${qc_dir} -o .
    """
}

/* -------------------------------------------------------------------------- */
workflow {

    file(params.outdir).mkdirs()

    /* --------------------------- Parameter checks -------------------------- */
    // For "real" runs, enforce these. (You can comment these out for toy demos.)
    require_param('samplesheet', params.samplesheet)
    require_param('ref', params.ref)
    require_param('gtf', params.gtf)
    require_param('transcripts', params.transcripts)

    /* ------------------------------- Reference ----------------------------- */
    ref_ch = Channel.fromPath(params.ref, checkIfExists: true)
    gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)
    tx_ch  = Channel.fromPath(params.transcripts, checkIfExists: true)

    star_index_ch   = ref_ch | STAR_INDEX
    salmon_index_ch = tx_ch  | SALMON_INDEX

    /* ------------------------------- Samples ------------------------------- */
    // Prefer samplesheet.csv: sample,read1,read2
    samples_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header:true)
        .map { row ->
            def sid = row.sample as String
            tuple(sid,
                  file(row.read1 as String),
                  file(row.read2 as String))
        }

    /* ------------------------------ fastp QC ------------------------------- */
    fastp_out_ch = samples_ch | FASTP_QC

    // Split trimmed reads back out
    trimmed_reads_ch = fastp_out_ch.map { sid, r1t, r2t, html, json -> tuple(sid, r1t, r2t) }

    /* --------------------------- STAR for tracks --------------------------- */
    aligned_ch = trimmed_reads_ch
        .combine(star_index_ch)
        .map { s_tup, idx -> tuple(s_tup[0], s_tup[1], s_tup[2], idx) }
        | STAR_ALIGN

    // Tracks
    aligned_ch | BAM_TO_BIGWIG

    /* ------------------------ featureCounts gene counts -------------------- */
    counts_ch = aligned_ch
        .combine(gtf_ch)
        .map { bam_tup, gtf -> tuple(bam_tup[0], bam_tup[1], gtf) }
        | FEATURECOUNTS_PER_SAMPLE

    counts_files_ch = counts_ch
        .map { sid, tsv, summary -> tsv }
        .collect()

    counts_files_ch | MERGE_COUNTS_MATRIX

    /* ------------------------------ Salmon quant --------------------------- */
    salmon_q_ch = trimmed_reads_ch
        .combine(salmon_index_ch)
        .map { s_tup, sidx -> tuple(s_tup[0], s_tup[1], s_tup[2], sidx) }
        | SALMON_QUANT

    // Stage quant.sf files as <sample>.quant.sf for merging
    salmon_q_ch
        .map { sid, qsf ->
            def staged = file("${sid}.quant.sf")
            staged.text = qsf.text
            staged
        }
        .collect()
        .set { staged_quant_files_ch }

    // Merge TPM matrix
    staged_quant_files_ch | MERGE_SALMON_TPM

    /* ------------------------------ MultiQC -------------------------------- */
    // Point MultiQC at the whole results directory (safe + simple)
    Channel.fromPath(params.outdir, type:'dir') | MULTIQC
}



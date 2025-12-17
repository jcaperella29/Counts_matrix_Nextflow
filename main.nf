nextflow.enable.dsl = 2

/*
  Portfolio-ready mini RNA-seq pipeline
  - QC: fastp
  - Genome alignment (for tracks + gene counts): STAR -> sorted BAM
  - Tracks: BigWig (deeptools bamCoverage)
  - Gene counts: wfeatureCounts (GTF)
  - Transcript quant: Salmon
  - Summary: MultiQC

  Required tools in container:
    fastp, STAR, samtools, subread(featureCounts), salmon, multiqc, deeptools

  Inputs: samplesheet.csv with columns: sample,read1,read2
*/

params.samplesheet = params.samplesheet ?: null
process MERGE_SALMON_TPM {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(quant_files)

    output:
    path "salmon_tpm_matrix.tsv"

    script:
    """
    python << 'EOF'
import os
import glob

qfiles = sorted(glob.glob('*.quant.sf'))
if not qfiles:
    raise SystemExit('No *.quant.sf files found for merging')

samples = []
rows = {}  # transcript -> {sample: TPM}

for q in qfiles:
    sample = os.path.basename(q).replace('.quant.sf', '')
    samples.append(sample)
    with open(q) as f:
        header = f.readline().strip().split('\\t')
        idx_name = header.index('Name')
        idx_tpm  = header.index('TPM')
        for line in f:
            parts = line.rstrip('\\n').split('\\t')
            tx = parts[idx_name]
            tpm = parts[idx_tpm]
            rows.setdefault(tx, {})[sample] = tpm

samples = sorted(samples)

with open('salmon_tpm_matrix.tsv', 'w') as out:
    out.write('transcript\\t' + '\\t'.join(samples) + '\\n')
    for tx in sorted(rows.keys()):
        out.write(tx)
        for s in samples:
            out.write('\\t' + rows[tx].get(s, '0'))
        out.write('\\n')
EOF
    """
}

params.input_dir   = params.input_dir   ?: 'data'
params.outdir      = params.outdir      ?: 'results'
params.threads     = params.threads     ?: 4

// Reference inputs
params.ref         = params.ref         ?: null
params.gtf         = params.gtf         ?: null
params.transcripts = params.transcripts ?: null

// Track settings (define defaults so BAM_TO_BIGWIG won't warn)
params.bw_binsize  = params.bw_binsize  ?: 10
params.bw_norm     = params.bw_norm     ?: 'CPM'

/* -------------------------------------------------------------------------- */
/* Helpers                                                                    */
/* -------------------------------------------------------------------------- */

def require_param(name, value) {
    if (!value) {
        throw new IllegalArgumentException("Missing required parameter --${name}")
    }
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
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}.flagstat.txt")

    script:
    """
    STAR \\
        --runThreadN ${params.threads} \\
        --genomeDir ${star_index_dir} \\
        --readFilesIn ${r1} ${r2} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${sample_id}.

# STAR output BAM (with your prefix) is:
#   ${sample_id}.Aligned.sortedByCoord.out.bam
       ls -lh ${sample_id}.* || true

  mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.bam
  samtools index ${sample_id}.bam
  samtools flagstat ${sample_id}.bam > ${sample_id}.flagstat.txt
 """

}

/* -------------------------------------------------------------------------- */
/* BAM_TO_BIGWIG – genome browser tracks via deeptools bamCoverage            */
/* -------------------------------------------------------------------------- */

process BAM_TO_BIGWIG {

  tag "${sample_id}"

  input:
    tuple val(sample_id), path(bam), path(bai)

  output:
    tuple val(sample_id), path("${sample_id}.bw")

  script:
  """
  mapped=\$(samtools view -c -F 4 ${bam})
  if [[ "\$mapped" -eq 0 ]]; then
    echo "No mapped reads for ${sample_id}; skipping bamCoverage"
    touch ${sample_id}.bw
  else
    bamCoverage \
      -b ${bam} \
      -o ${sample_id}.bw \
      --binSize ${params.bw_binsize} \
      --normalizeUsing ${params.bw_norm} \
      -p ${params.threads}
  fi
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
    tuple val(sample_id),
          path("${sample_id}.counts.tsv"),
          path("${sample_id}.featureCounts.summary")

    script:
    """
    featureCounts \\
      -T ${params.threads} \\
      -a ${gtf} \\
      -o ${sample_id}.raw_counts.txt \\
      ${bam}

    cp ${sample_id}.raw_counts.txt.summary ${sample_id}.featureCounts.summary

    awk 'BEGIN{OFS="\\t"} \
         /^#/ {next} \
         NR==2 {next} \
         NR>2 {print \$1, \$7}' \
         ${sample_id}.raw_counts.txt \
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
  python - <<'PY'
  import glob, os, csv

  files = sorted(glob.glob("*.counts.tsv"))
  if not files:
      raise SystemExit("No count files found.")

  all_counts = {}
  features = set()

  for fname in files:
      sid = os.path.basename(fname).split(".")[0]
      all_counts[sid] = {}
      with open(fname, newline="") as f:
          for row in csv.reader(f, delimiter="\\t"):
              if not row:
                  continue
              feat = row[0]
              count = row[-1]
              all_counts[sid][feat] = count
              features.add(feat)

  sids = sorted(all_counts.keys())
  features = sorted(features)

  with open("counts_matrix.tsv", "w", newline="") as out:
      out.write("feature\\t" + "\\t".join(sids) + "\\n")
      for feat in features:
          out.write(feat + "\\t" + "\\t".join(all_counts[sid].get(feat, "0") for sid in sids) + "\\n")
  PY
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
    tuple val(sample_id), path("${sample_id}.quant.sf")

  script:
    """
    salmon quant \
      -i ${salmon_index} \
      -l A \
      -1 ${r1} -2 ${r2} \
      -p ${params.threads} \
      --validateMappings \
      -o .

    cp quant.sf ${sample_id}.quant.sf
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

    require_param('samplesheet', params.samplesheet)
    require_param('ref', params.ref)
    require_param('gtf', params.gtf)
    require_param('transcripts', params.transcripts)

    ref_ch = Channel.fromPath(params.ref, checkIfExists: true)
    gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)
    tx_ch  = Channel.fromPath(params.transcripts, checkIfExists: true)

    star_index_ch   = ref_ch | STAR_INDEX
    salmon_index_ch = tx_ch  | SALMON_INDEX

    samples_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header:true)
        .map { row ->
            def sid = row.sample as String
            tuple(sid, file(row.read1 as String), file(row.read2 as String))
        }

    fastp_out_ch = samples_ch | FASTP_QC

    trimmed_reads_ch = fastp_out_ch.map { sid, r1t, r2t, html, json -> tuple(sid, r1t, r2t) }

    /*
      KEY FIX: After `.combine()`, Nextflow will pass tuple elements as separate
      args into the closure. So the map closure MUST accept 4 args here.
    */

   aligned_ch = trimmed_reads_ch
  .combine(star_index_ch)
  | STAR_ALIGN



    aligned_ch
       .map { sid, bam, bai, flagstat -> tuple(sid, bam, bai) }
       | BAM_TO_BIGWIG


    counts_ch = aligned_ch
  .combine(gtf_ch)
  .map { sid, bam, bai, flagstat, gtf ->
      tuple(sid, bam, gtf)
  }
  | FEATURECOUNTS_PER_SAMPLE



    counts_files_ch = counts_ch
        .map { sid, tsv, summary -> tsv }
        .collect()

    counts_files_ch | MERGE_COUNTS_MATRIX

    
   salmon_in_ch = trimmed_reads_ch
  .combine(salmon_index_ch)

salmon_q_ch = salmon_in_ch | SALMON_QUANT





   


    salmon_quant_files_ch = salmon_q_ch
        .map { sid, qsf -> qsf }
        .collect()

    salmon_quant_files_ch | MERGE_SALMON_TPM

    Channel.fromPath(params.outdir, type:'dir') | MULTIQC
}






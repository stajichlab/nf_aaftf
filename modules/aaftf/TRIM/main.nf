process AAFTF_TRIM {
    tag   { sample }
    label 'aaftf_trim'
    publishDir "${params.outdir}/trim", mode: 'copy', pattern: '_1P.fastq.gz|_2P.fastq.gz|_MG.fastq.gz'

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    tuple val(sample), path("${sample}_1P.fastq.gz"), path("${sample}_2P.fastq.gz"), path("${sample}_MG.fastq.gz"), emit: trimmed
    path "trim_*.json", emit: json, optional: true

    script:
    // Pass-1 trimer: fastp (default, adds dedup+merge+3' qual trim) or
    // trimmomatic (separate adaptor/leading/trailing/sliding-window knobs).
    // Only fastp emits a merged read set (_MG); trimmomatic yields an empty one.
    def m1      = params.trim_method1 ?: 'fastp'
    def fastp   = m1.toLowerCase() == 'fastp'
    def dedup   = fastp && params.dedup ? '--dedup' : ''
    def merge   = fastp ? '--merge' : ''
    def cuttail = (fastp && params.PHRED > 0) ? '--cuttail' : ''
    def trimmomatic_args = fastp ? '' : """
        --trimmomatic_adaptors ${params.trimmomatic_adaptors} \\
        --trimmomatic_clip ${params.trimmomatic_clip} \\
        --trimmomatic_leadingwindow ${params.trimmomatic_leadingwindow} \\
        --trimmomatic_trailingwindow ${params.trimmomatic_trailingwindow} \\
        --trimmomatic_slidingwindow ${params.trimmomatic_slidingwindow} \\
        --trimmomatic_quality ${params.trimmomatic_quality}
        """
    def merged_out = fastp ? "mv ${sample}_pass1_MG.fastq.gz ${sample}_MG.fastq.gz" :
                             "gzip -c /dev/null > ${sample}_MG.fastq.gz"
    """
    # Pass 1: ${m1} — QC + adaptor trimming. fastp adds dedup + merge + 3' qual
    # (cuttail using PHRED); trimmomatic uses its own leading/trailing/sliding
    # window settings.
    AAFTF trim --method ${m1} \
        ${dedup} ${merge} ${cuttail} --minlen ${params.minlen} \
        -c ${task.cpus} -m 8 ${trimmomatic_args} \
        --left ${reads1} --right ${reads2} \
        -o ${sample}_pass1

    # Pass 2: bbduk — a second adapter/quality pass (JGI-style BBDuk QC).
    # AAFTF would run `AAFTF trim --method bbduk` on the paired pass1 files, but
    # bbduk 39.80 in the image crashes its PairStreamer on variable-length paired
    # reads ("List size mismatch"), dropping all reads. Feed it an INTERLEAVED
    # single file instead (same effective result), then de-interleave.
    shuffle.sh in1=${sample}_pass1_1P.fastq.gz in2=${sample}_pass1_2P.fastq.gz out=${sample}_pass1_ivl.fq
    bbduk.sh -Xmx${params.trim_memory}g \\
        ref=adapters t=${task.cpus} ktrim=r k=23 mink=11 \\
        minlen=${params.minlen} hdist=1 maq=10 ftm=5 tpe tbo \\
        overwrite=true in=${sample}_pass1_ivl.fq interleaved=true \\
        out=${sample}_ivl.fq
    reformat.sh in=${sample}_ivl.fq out1=${sample}_1P.fastq.gz out2=${sample}_2P.fastq.gz

    # Merge-only reads (from fastp pass 1) pass forward as the single-end set.
    # Trimmomatic has no merge step, so an empty _MG is emitted for it.
    ${merged_out}
    """

    stub:
    """
    touch ${sample}_1P.fastq.gz
    touch ${sample}_2P.fastq.gz
    touch ${sample}_MG.fastq.gz
    touch trim_${sample}.json
    """
}

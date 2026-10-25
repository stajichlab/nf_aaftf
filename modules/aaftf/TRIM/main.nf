process AAFTF_TRIM {
    tag   { sample }
    label 'aaftf_trim'
    // Previous pattern ('_1P.fastq.gz|_2P.fastq.gz|_MG.fastq.gz') is not a
    // valid Nextflow glob (no leading wildcard, '|' isn't glob alternation),
    // so nothing was ever actually copied to outdir/trim -- confirmed empty
    // 2026-09-26. publishDir with no pattern copies exactly the declared
    // `output` paths, which is already just these three files.
    publishDir "${params.outdir}/trim", mode: 'copy'

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
    // BBTools picks its own Java heap when none is given (4.5 GB inside a 32 GB job). The interleave step used
    // shuffle.sh, which holds every read in memory: it ran out of heap on large samples even at 24 GB, and
    // the Java crash left the shell running. reformat.sh interleaves the pairs as a stream, so it needs almost
    // no memory. Give the BBTools steps 75% of the task memory anyway.
    // The interleaved intermediates are gzipped (BBTools picks gzip from the .gz extension): plain, one was
    // 15.5 GB for a pair of 1.5 and 1.8 GB read files, about 69% of the work dir.
    def heap = Math.max(4, (task.memory.toGiga() * 0.75) as int)
    """
    set -euo pipefail
    # (the aaftf profile sets shell=/bin/bash without -ue; without this a crashed step still exits 0)
    trap 'echo "AAFTF_TRIM: ERROR: command failed at line \$LINENO of .command.sh (exit \$?)" >&2' ERR

    # Pass 1: ${m1} — QC + adaptor trimming. fastp adds dedup + merge + 3' qual
    # (cuttail using PHRED); trimmomatic uses its own leading/trailing/sliding
    # window settings.
    AAFTF trim --method ${m1} \
        ${dedup} ${merge} ${cuttail} --minlen ${params.minlen} \
        -c ${task.cpus} -m 8 ${trimmomatic_args} \
        -1 ${reads1} -2 ${reads2} \
        -o ${sample}_pass1

    # Pass 2: bbduk — a second adapter/quality pass (JGI-style BBDuk QC).
    # AAFTF would run `AAFTF trim --method bbduk` on the paired pass1 files, but
    # bbduk 39.80 in the image crashes its PairStreamer on variable-length paired
    # reads ("List size mismatch"), dropping all reads. Feed it an INTERLEAVED
    # single file instead (same effective result), then de-interleave.
    reformat.sh -Xmx${heap}g in1=${sample}_pass1_1P.fastq.gz in2=${sample}_pass1_2P.fastq.gz out=${sample}_pass1_ivl.fq.gz
    # BBTools can return 0 after a Java crash, so check that each step really wrote its output.
    test -s ${sample}_pass1_ivl.fq.gz
    bbduk.sh -Xmx${params.trim_memory}g \\
        ref=adapters t=${task.cpus} ktrim=r k=23 mink=11 \\
        minlen=${params.minlen} hdist=1 maq=10 ftm=5 tpe tbo \\
        overwrite=true in=${sample}_pass1_ivl.fq.gz interleaved=true \\
        out=${sample}_ivl.fq.gz
    test -s ${sample}_ivl.fq.gz
    reformat.sh -Xmx${heap}g in=${sample}_ivl.fq.gz out1=${sample}_1P.fastq.gz out2=${sample}_2P.fastq.gz
    test -s ${sample}_1P.fastq.gz && test -s ${sample}_2P.fastq.gz

    # Merge-only reads (from fastp pass 1) pass forward as the single-end set.
    # Trimmomatic has no merge step, so an empty _MG is emitted for it.
    ${merged_out}
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'fastp=fastp --version' \\
        'bbduk=bbduk.sh --version 2>&1 | grep -i '\\''BBTools version'\\''' \\
        'pigz=pigz --version'
    """

    stub:
    """
    touch ${sample}_1P.fastq.gz
    touch ${sample}_2P.fastq.gz
    touch ${sample}_MG.fastq.gz
    touch trim_${sample}.json
    """
}

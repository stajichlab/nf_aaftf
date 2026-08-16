process AFFTF_TRIM {
    tag   { sample }
    label 'aaftf_trim'
    publishDir "${params.outdir}/trim", mode: 'copy', pattern: '_1P.fastq.gz|_2P.fastq.gz|_MG.fastq.gz'

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    tuple val(sample), path("${sample}_1P.fastq.gz"), path("${sample}_2P.fastq.gz"), path("${sample}_MG.fastq.gz"), emit: trimmed
    path "trim_*.json", emit: json, optional: true

    script:
    def dedup   = params.dedup ? '--dedup' : ''
    def cuttail = params.PHRED > 0 ? '--cuttail' : ''
    """
    # Pass 1: fastp — dedup + merge + 3' quality trimming (cuttail using PHRED).
    AAFTF trim --method fastp \
        ${dedup} --merge ${cuttail} --minlen ${params.minlen} \
        -c ${task.cpus} -m 8 \
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

    # Merge-only reads (from pass 1) pass forward as the single-end set.
    mv ${sample}_pass1_MG.fastq.gz ${sample}_MG.fastq.gz
    """

    stub:
    """
    touch ${sample}_1P.fastq.gz
    touch ${sample}_2P.fastq.gz
    touch ${sample}_MG.fastq.gz
    touch trim_${sample}.json
    """
}

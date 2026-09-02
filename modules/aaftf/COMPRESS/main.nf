process COMPRESS {
    tag   { sample }
    label 'aaftf_lite'
    publishDir "${params.outdir}/sort", mode: 'copy', pattern: '*.sorted.fasta.gz'

    input:
    tuple val(sample), path(asm)

    output:
    tuple val(sample), path("${sample}.sorted.fasta.gz"), emit: compressed

    script:
    """
    # bgzip (htslib), not plain gzip: keeps the final assembly block-compressed
    # and indexable (samtools faidx / bcftools) while remaining gzip-compatible.
    bgzip -c -@ ${task.cpus} ${asm} > ${sample}.sorted.fasta.gz
    """

    stub:
    """
    gzip -c ${asm} > ${sample}.sorted.fasta.gz
    """
}

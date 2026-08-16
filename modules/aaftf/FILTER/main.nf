process FILTER {
    tag   { sample }
    label 'aaftf_trim'
    publishDir "${params.outdir}/filter", mode: 'copy', pattern: '_filtered_1.fastq.gz|_filtered_2.fastq.gz|_filtered_U.fastq.gz'

    input:
    tuple val(sample), path(reads1), path(reads2), path(merged)

    output:
    tuple val(sample), path("${sample}_filtered_1.fastq.gz"), path("${sample}_filtered_2.fastq.gz"), path("${sample}_filtered_U.fastq.gz"), emit: filtered

    script:
    """
    # Paired-end reads: filter against phiX + common contaminant/vector DBs.
    AAFTF filter --aligner ${params.filter_aligner} \
        -c ${task.cpus} -m 16 --AAFTF_DB /opt/aaaftf_db \
        --left ${reads1} --right ${reads2} \
        -o ${sample}

    # Merged (single-end) reads: same filter run as single-end so merged reads
    # survive too; AAFTF emits them as _filtered_U.
    AAFTF filter --aligner ${params.filter_aligner} \
        -c ${task.cpus} -m 16 --AAFTF_DB /opt/aaaftf_db \
        --left ${merged} \
        -o ${sample}
    """

    stub:
    """
    touch ${sample}_filtered_1.fastq.gz
    touch ${sample}_filtered_2.fastq.gz
    touch ${sample}_filtered_U.fastq.gz
    """
}

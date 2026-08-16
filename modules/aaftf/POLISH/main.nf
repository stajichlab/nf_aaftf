process POLISH {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/polish", mode: 'copy', pattern: '*.polished.fasta'

    input:
    tuple val(sample), path(asm), path(filtered_1), path(filtered_2)

    output:
    tuple val(sample), path("${sample}.polished.fasta"), emit: polished

    script:
    """
    # Polish with short reads. POLCA is the default in the source framework and
    # is robust to indel error; map the ORIGINAL filtered reads back to the
    # assembly (post-rmdup) for error correction.
    AAFTF polish --method ${params.polisher} \
        -i ${asm} -o ${sample}.polished.fasta \
        -c ${task.cpus} -m ${params.spades_memory} \
        --left ${filtered_1} --right ${filtered_2} \
        -w polish_work_${sample}
    """

    stub:
    """
    cp ${asm} ${sample}.polished.fasta
    """
}

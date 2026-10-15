process POLISH {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/polish", mode: 'copy', pattern: '*.polished.fasta.gz'

    input:
    tuple val(sample), path(asm), path(filtered_1), path(filtered_2)

    output:
    tuple val(sample), path("${sample}.polished.fasta"), emit: polished
    path "${sample}.polished.fasta.gz"

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
    gzip -c ${sample}.polished.fasta > ${sample}.polished.fasta.gz
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version' \\
        'pilon=pilon --version 2>&1 | head -1' \\
        'bwa=bwa 2>&1 | grep -i '\\''^Version'\\''' \\
        'samtools=samtools --version | head -1' \\
        "polisher=echo ${params.polisher}"
    """

    stub:
    """
    cp ${asm} ${sample}.polished.fasta
    gzip -c ${sample}.polished.fasta > ${sample}.polished.fasta.gz
    """
}

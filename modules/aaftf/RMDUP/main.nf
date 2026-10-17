process RMDUP {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/rmdup", mode: 'copy', pattern: '*.rmdup.fasta.gz'

    input:
    tuple val(sample), path(asm)

    output:
    tuple val(sample), path("${sample}.rmdup.fasta"), emit: rmdup
    path "${sample}.rmdup.fasta.gz"

    script:
    """
    # Remove redundant (duplicate) contigs. min length keeps only contigs
    # >= params.min_contig_len (post-assembly contig cleanup).
    AAFTF rmdup -i ${asm} -o ${sample}.rmdup.fasta \
        -c ${task.cpus} -ml ${params.min_contig_len}
    gzip -c ${sample}.rmdup.fasta > ${sample}.rmdup.fasta.gz
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'minimap2=minimap2 --version' \\
        'rmdup_tools=echo minimap2'
    """

    stub:
    """
    cp ${asm} ${sample}.rmdup.fasta
    gzip -c ${sample}.rmdup.fasta > ${sample}.rmdup.fasta.gz
    """
}

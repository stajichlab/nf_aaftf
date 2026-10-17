process SORT {
    tag   { sample }
    label 'aaftf_lite'
    // Final assembly is published bgzip-compressed by COMPRESS (see main.nf);
    // this uncompressed copy stays in work/ only, for ASSESS/DEPTH to consume.

    input:
    tuple val(sample), path(asm)

    output:
    tuple val(sample), path("${sample}.sorted.fasta"), emit: sorted

    script:
    """
    # Sort contigs by length (largest first), rename headers to the sample name,
    # and drop contigs shorter than params.min_contig_len.
    AAFTF sort -i ${asm} -o ${sample}.sorted.fasta \
        -ml ${params.min_contig_len} -n ${sample}
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1'
    """

    stub:
    """
    cp ${asm} ${sample}.sorted.fasta
    """
}

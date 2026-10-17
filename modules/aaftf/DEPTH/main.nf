process DEPTH {
    tag   { sample }
    label 'aaftf_lite'
    publishDir "${params.outdir}/depth", mode: 'copy', pattern: '*.depth.txt'

    input:
    tuple val(sample), path(asm), path(filtered_1), path(filtered_2)

    output:
    tuple val(sample), path("${sample}.depth.txt"), emit: depth

    script:
    """
    # Map reads back with minimap2 and compute per-contig depth (flags outlier
    # contigs as possible contaminants / organellar). Disable plots to avoid the
    # matplotlib dependency in headless runs.
    AAFTF depth -i ${asm} -o ${sample}.depth.txt \
        -1 ${filtered_1} -2 ${filtered_2} --no-plot \
        -c ${task.cpus} -w depth_work_${sample}
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'minimap2=minimap2 --version' \\
        'samtools=samtools --version | head -1' \\
        'mosdepth=mosdepth --version'
    """

    stub:
    """
    cat > ${sample}.depth.txt <<'EOF'
contig	depth
${sample}_contig_1	42.0
EOF
    """
}

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
        -l ${filtered_1} -r ${filtered_2} --no-plot \
        -c ${task.cpus} -w depth_work_${sample}
    """

    stub:
    """
    cat > ${sample}.depth.txt <<'EOF'
contig	depth
${sample}_contig_1	42.0
EOF
    """
}

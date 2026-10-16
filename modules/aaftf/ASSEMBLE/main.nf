process ASSEMBLE {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/asm", mode: 'copy', pattern: "*.${params.assembler}.fasta.gz"

    input:
    tuple val(sample), path(filtered_1), path(filtered_2), path(filtered_U)

    output:
    tuple val(sample), path("${sample}.${params.assembler}.fasta"), emit: assembly
    path "${sample}.${params.assembler}.fasta.gz", emit: assembly_gz

    script:
    // Assembler switch: spades (default) | megahit | unicycler. (AAFTF >= 0.7 dropped dipspades.)
    // AAFTF adds --careful + --cov-cutoff auto for spades; megahit/unicycler get
    // their own defaults. Extra args pass through via params.assembler_args.
    // `--merged` gets the filtered merged reads as the single-end input.
    def extra = params.assembler_args ? "--assembler_args '${params.assembler_args}'" : ''
    """
    AAFTF assemble --method ${params.assembler} \
        -c ${task.cpus} -m ${params.spades_memory} \
        -1 ${filtered_1} -2 ${filtered_2} \
        --merged ${filtered_U} \
        -o ${sample}.${params.assembler}.fasta \
        -w spades_work_${sample} ${extra}
    gzip -c ${sample}.${params.assembler}.fasta > ${sample}.${params.assembler}.fasta.gz
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'spades=spades.py --version' \\
        'megahit=megahit --version' \\
        "assembler=echo ${params.assembler}"
    """

    stub:
    """
    cat > ${sample}.${params.assembler}.fasta <<'FASTA'
>${sample}_contig_1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
FASTA
    gzip -c ${sample}.${params.assembler}.fasta > ${sample}.${params.assembler}.fasta.gz
    """
}

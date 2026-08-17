process ASSEMBLE {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/asm", mode: 'copy', pattern: "*.${params.assembler}.fasta"

    input:
    tuple val(sample), path(filtered_1), path(filtered_2), path(filtered_U)

    output:
    tuple val(sample), path("${sample}.${params.assembler}.fasta"), emit: assembly
    path "${sample}.${params.assembler}.fasta.gz", emit: assembly_gz, optional: true

    script:
    // Assembler switch: spades (default) | dipspades | megahit | unicycler.
    // AAFTF adds --careful + --cov-cutoff auto for spades; megahit/unicycler get
    // their own defaults. Extra args pass through via params.assembler_args.
    // `--merged` gets the filtered merged reads as the single-end input.
    def extra = params.assembler_args ? "--assembler_args '${params.assembler_args}'" : ''
    """
    AAFTF assemble --method ${params.assembler} \
        -c ${task.cpus} -m ${params.spades_memory} \
        --left ${filtered_1} --right ${filtered_2} \
        --merged ${filtered_U} \
        -o ${sample}.${params.assembler}.fasta \
        -w spades_work_${sample} ${extra}
    """

    stub:
    """
    cat > ${sample}.${params.assembler}.fasta <<'FASTA'
>${sample}_contig_1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
FASTA
    """
}

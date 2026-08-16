process ASSEMBLE {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/asm", mode: 'copy', pattern: '*.spades.fasta'

    input:
    tuple val(sample), path(filtered_1), path(filtered_2), path(filtered_U)

    output:
    tuple val(sample), path("${sample}.spades.fasta"), emit: assembly
    path "${sample}.spades.fasta.gz", emit: assembly_gz, optional: true

    script:
    // SPAdes (short-read Illumina). AAFTF adds --careful + --cov-cutoff auto.
    // `--merged` gets the filtered merged reads as the single-end input.
    """
    AAFTF assemble --method ${params.assembler} \
        -c ${task.cpus} -m ${params.spades_memory} \
        --left ${filtered_1} --right ${filtered_2} \
        --merged ${filtered_U} \
        -o ${sample}.spades.fasta \
        -w spades_work_${sample}
    """

    stub:
    """
    cat > ${sample}.spades.fasta <<'FASTA'
>${sample}_contig_1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
FASTA
    """
}

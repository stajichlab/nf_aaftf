process VECSCREEN {
    tag   { sample }
    label 'aaftf_lite'
    publishDir "${params.outdir}/vecscreen", mode: 'copy', pattern: '*.vecscreen.fasta'

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("${sample}.vecscreen.fasta"), emit: vecscreen

    script:
    """
    # BLASTN screen of contigs against UniVec + contaminant DBs (vecscreen mode).
    # Note: this runs fully inside the AAFTF container — no external FCS tool.
    AAFTF vecscreen -c ${task.cpus} \\
        --AAFTF_DB /opt/aaftf_db \\
        -i ${assembly} -o ${sample}.vecscreen.fasta
    """

    stub:
    """
    cp ${assembly} ${sample}.vecscreen.fasta
    """
}

process VECSCREEN {
    tag   { sample }
    label 'aaftf_lite'
    publishDir "${params.outdir}/vecscreen", mode: 'copy', pattern: '*.vecscreen.fasta.gz'

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("${sample}.vecscreen.fasta"), emit: vecscreen
    path "${sample}.vecscreen.fasta.gz"

    script:
    """
    # BLASTN screen of contigs against UniVec + contaminant DBs (vecscreen mode).
    # Note: this runs fully inside the AAFTF container — no external FCS tool.
    AAFTF vecscreen -c ${task.cpus} \\
        --AAFTF_DB /opt/aaftf_db \\
        -i ${assembly} -o ${sample}.vecscreen.fasta
    gzip -c ${sample}.vecscreen.fasta > ${sample}.vecscreen.fasta.gz
    """

    stub:
    """
    cp ${assembly} ${sample}.vecscreen.fasta
    gzip -c ${sample}.vecscreen.fasta > ${sample}.vecscreen.fasta.gz
    """
}

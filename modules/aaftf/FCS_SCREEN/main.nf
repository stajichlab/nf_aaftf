process FCS_SCREEN {
    tag   { sample }
    label 'aaftf_lite'
    publishDir "${params.outdir}/fcs_screen", mode: 'copy', pattern: '*.fcs_screen.fasta'

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("${sample}.fcs_screen.fasta"), emit: screened

    script:
    // NCBI FCS adaptor screening — an alternative to vecscreen for vector
    // and common contamination screening. Uses --euk (default) or --prok.
    def mode = (params.getOrDefault('fcs_screen_prok', false) as String).toBoolean() ? '--prok' : '--euk'
    """
    AAFTF fcs_screen -i ${assembly} -o ${sample}.fcs_screen.fasta ${mode}
    """

    stub:
    """
    cp ${assembly} ${sample}.fcs_screen.fasta
    """
}

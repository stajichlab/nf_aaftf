process FCS_SCREEN {
    tag   { sample }
    label 'aaftf_native'
    publishDir "${params.outdir}/fcs_screen", mode: 'copy', pattern: '*.fcs_screen.fasta'
    publishDir "${params.outdir}/fcs_screen", mode: 'copy', pattern: '*.fcs_adaptor_report.txt'

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("${sample}.fcs_screen.fasta"), emit: screened
    path "${sample}.fcs_adaptor_report.txt"

    script:
    // NCBI FCS adaptor screening — an alternative to vecscreen for vector
    // and common contamination screening. Uses --euk (default) or --prok.
    //
    // AAFTF fcs_screen just shells out to NCBI's run_fcsadaptor.sh, which
    // itself launches a SECOND singularity/apptainer container running
    // /app/fcs/bin/av_screen_x inside fcs-adaptor.sif. Nested containers
    // don't work from inside our own AAFTF SIF (see withLabel 'aaftf_native',
    // container = false), and run_fcsadaptor.sh's own singularity branch is
    // broken anyway (it does `singularity run $CONTAINER --bind ...` with
    // $CONTAINER unset — a docker-branch variable leftover). So skip both
    // AAFTF's wrapper and NCBI's wrapper script and call av_screen_x inside
    // the cached fcs-adaptor SIF directly with apptainer exec, replicating
    // exactly what run_fcsadaptor.sh's (working) docker branch does. This
    // also sidesteps an AAFTF fcs_screen.py bug where a trailing
    // `os.mkdir(cleaned_sequences)` crashes because av_screen_x already
    // created that directory, even though screening itself succeeded.
    def mode = (params.getOrDefault('fcs_screen_prok', false) as String).toBoolean() ? '--prok' : '--euk'
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load apptainer

    EXPANDED_FASTA=\$(readlink -f ${assembly})
    FASTA_DIR=\$(dirname "\$EXPANDED_FASTA")
    FASTA_NAME=\$(basename "\$EXPANDED_FASTA")
    OUTDIR=\$(readlink -f fcs_screen_out)
    mkdir -p "\$OUTDIR"

    apptainer exec --bind "\$FASTA_DIR":/sample-volume --bind "\$OUTDIR":/output-volume \\
        ${params.fcs_adaptor_sif} \\
        /app/fcs/bin/av_screen_x -o /output-volume ${mode} /sample-volume/\$FASTA_NAME

    cp "\$OUTDIR/cleaned_sequences/\$FASTA_NAME" ${sample}.fcs_screen.fasta
    cp "\$OUTDIR/fcs_adaptor_report.txt" ${sample}.fcs_adaptor_report.txt
    """

    stub:
    """
    cp ${assembly} ${sample}.fcs_screen.fasta
    touch ${sample}.fcs_adaptor_report.txt
    """
}

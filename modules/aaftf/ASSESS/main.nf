process ASSESS {
    tag   { sample }
    label 'aaftf_lite'
    publishDir "${params.outdir}/assess", mode: 'copy', pattern: '*.stats.txt'

    input:
    tuple val(sample), path(asm)

    output:
    tuple val(sample), path("${sample}.stats.txt"), emit: stats

    script:
    """
    # Assembly statistics (N50, contig count, GC%, etc.) written to a report.
    AAFTF assess -i ${asm} -r ${sample}.stats.txt
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version'
    """

    stub:
    """
    cat > ${sample}.stats.txt <<'EOF'
CONTIG COUNT  =  10
TOTAL LENGTH  =  1000000
         N50  =  150000
         L50  =  3
           GC %  =  50.0
EOF
    """
}

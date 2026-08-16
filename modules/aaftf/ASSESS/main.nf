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

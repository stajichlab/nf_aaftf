// Merge several sequencing libraries of ONE sample into a single read pair.
// The sample sheet may list a sample on several rows (one row per library).
// gzip members can be concatenated, so `cat` joins them without recompression.
// The component order is the order given by main.nf (sorted by read_1 file name).
process MERGE_LIBRARIES {
    tag   { sample }
    label 'aaftf_lite'
    publishDir "${params.outdir}/libraries", mode: 'copy', pattern: '*.libraries.tsv'

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    tuple val(sample), path("${sample}_R1.merged.fq.gz"), path("${sample}_R2.merged.fq.gz"), emit: merged
    path "${sample}.libraries.tsv", emit: manifest

    script:
    """
    set -euo pipefail
    cat ${reads1} > ${sample}_R1.merged.fq.gz
    cat ${reads2} > ${sample}_R2.merged.fq.gz
    gzip -t ${sample}_R1.merged.fq.gz
    gzip -t ${sample}_R2.merged.fq.gz
    # Mates must stay in step: equal line counts for R1 and R2.
    n1=\$(gzip -dc ${sample}_R1.merged.fq.gz | wc -l)
    n2=\$(gzip -dc ${sample}_R2.merged.fq.gz | wc -l)
    if [ "\$n1" != "\$n2" ]; then
        echo "ERROR: ${sample}: merged R1 has \$n1 lines but R2 has \$n2" >&2
        exit 1
    fi
    printf 'sample\\tread_1\\tread_2\\n' > ${sample}.libraries.tsv
    paste <(printf '%s\\n' ${reads1}) <(printf '%s\\n' ${reads2}) | awk -v s=${sample} 'BEGIN{OFS="\\t"}{print s,\$1,\$2}' >> ${sample}.libraries.tsv
    """

    stub:
    """
    printf '' | gzip -c > ${sample}_R1.merged.fq.gz
    printf '' | gzip -c > ${sample}_R2.merged.fq.gz
    printf 'sample\\tread_1\\tread_2\\n' > ${sample}.libraries.tsv
    """
}

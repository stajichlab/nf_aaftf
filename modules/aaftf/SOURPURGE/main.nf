process SOURPURGE {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/sourpurge", mode: 'copy', pattern: '*.sourpurge.fasta.gz'

    input:
    tuple val(sample), path(asm), val(taxonid)

    output:
    tuple val(sample), path("${sample}.sourpurge.fasta"), emit: clean
    path "${sample}.sourpurge.fasta.gz"

    script:
    // Sourmash-based contamination purge — an alternative or complement to
    // FCS-GX. Compares contig k-mers against a sourmash LCA taxonomy DB and
    // purges contigs whose taxonomy doesn't match the expected phylum.
    // Resolve phylum NAME from the organism taxid via taxonkit.
    """
    TAXONKIT_DB=${params.taxondb}
    PHYLUM=\$(echo "${taxonid}" | taxonkit --data-dir \$TAXONKIT_DB lineage \\
        | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result \\
        | cut -f3)
    if [ -z "\$PHYLUM" ] || [ "\$PHYLUM" = "No rank" ]; then
        echo "SOURPURGE: could not resolve phylum name from taxid ${taxonid}" >&2
        exit 1
    fi
    echo "SOURPURGE: ${sample} taxonid=${taxonid} -> phylum='\$PHYLUM'"

    # AAFTF >= 0.7: `sourpurge` was split into `taxpurge` (sourmash taxonomy screen) and
    # `covpurge` (bwa low-coverage screen, needs reads). This step is the taxonomy screen.
    # Needs the sourmash databases in the AAFTF_DB folder (`AAFTF database taxpurge`).
    export AAFTF_DB=/opt/aaftf_db
    AAFTF taxpurge \\
        -i ${asm} -o ${sample}.sourpurge.fasta \\
        -p "\$PHYLUM" -c ${task.cpus}
    gzip -c ${sample}.sourpurge.fasta > ${sample}.sourpurge.fasta.gz
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'sourmash=sourmash --version' \\
        'taxonkit=taxonkit version' \\
        "aaftf_db=readlink -f ${params.aaftf_db}"
    """

    stub:
    """
    cp ${asm} ${sample}.sourpurge.fasta
    gzip -c ${sample}.sourpurge.fasta > ${sample}.sourpurge.fasta.gz
    """
}

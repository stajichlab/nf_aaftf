process SOURPURGE {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/sourpurge", mode: 'copy', pattern: '*.sourpurge.fasta'

    input:
    tuple val(sample), path(asm), val(taxonid)

    output:
    tuple val(sample), path("${sample}.sourpurge.fasta"), emit: clean

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

    AAFTF sourpurge \\
        -i ${asm} -o ${sample}.sourpurge.fasta \\
        -p "\$PHYLUM" -c ${task.cpus} \\
        --AAFTF_DB /opt/aaftf_db
    """

    stub:
    """
    cp ${asm} ${sample}.sourpurge.fasta
    """
}

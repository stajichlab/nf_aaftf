process SOURPURGE {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/sourpurge", mode: 'copy', pattern: '*.sourpurge.fasta.gz'

    input:
    tuple val(sample), path(asm), val(phylum_name)

    output:
    tuple val(sample), path("${sample}.sourpurge.fasta"), emit: clean
    path "${sample}.sourpurge.fasta.gz"

    script:
    // Sourmash-based contamination purge — an alternative or complement to
    // FCS-GX. Compares contig k-mers against a sourmash LCA taxonomy DB and
    // purges contigs whose taxonomy doesn't match the expected phylum.
    // The phylum NAME comes from TAXON_PHYLUM (taxonkit runs in its own container).
    """
    PHYLUM="${phylum_name}"
    if [ -z "\$PHYLUM" ]; then
        echo "SOURPURGE: no phylum name resolved for ${sample}" >&2
        exit 1
    fi
    echo "SOURPURGE: ${sample} phylum='\$PHYLUM'"

    # AAFTF >= 0.7: `sourpurge` was split into `taxpurge` (sourmash taxonomy screen) and
    # `covpurge` (bwa low-coverage screen, needs reads). This step is the taxonomy screen.
    # Needs the sourmash databases in the AAFTF_DB folder (`AAFTF database taxpurge`).
    export AAFTF_DB=${params.getOrDefault('aaftf_db_shared', null) ? '/opt/aaftf_db:/opt/aaftf_db_shared' : '/opt/aaftf_db'}
    AAFTF taxpurge \\
        -i ${asm} -o ${sample}.sourpurge.fasta \\
        -p "\$PHYLUM" -c ${task.cpus}
    gzip -c ${sample}.sourpurge.fasta > ${sample}.sourpurge.fasta.gz
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'sourmash=sourmash --version' \\
        "aaftf_db=readlink -f ${params.aaftf_db}"
    """

    stub:
    """
    cp ${asm} ${sample}.sourpurge.fasta
    gzip -c ${sample}.sourpurge.fasta > ${sample}.sourpurge.fasta.gz
    """
}

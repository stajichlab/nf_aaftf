process CONTAM_CLEAN {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/contam_clean", mode: 'copy', pattern: '*.contam_clean.fasta.gz'

    input:
    tuple val(sample), path(asm), val(phylum)

    output:
    tuple val(sample), path("${sample}.contam_clean.fasta"), emit: clean
    path "${sample}.contam_clean.fasta.gz"

    script:
    """
    # Optional NCBI FCS-GX contamination purge. `phylum` is the numeric PHYLUM
    # taxid FCS-GX -t expects. TAXON_PHYLUM resolves it with taxonkit in its own
    # container (AAFTF >= 0.7 no longer ships taxonkit).
    if ! [[ "${phylum}" =~ ^[0-9]+\$ ]]; then
        echo "CONTAM_CLEAN: ERROR: phylum '${phylum}' is not numeric; cannot run fcs_gx_purge for ${sample}" >&2
        exit 1
    fi
    echo "CONTAM_CLEAN: ${sample} phylum_taxid=${phylum}"

    # The FCS-GX database is large; stage it to node-local scratch once per sample.
    STAGE="\${SCRATCH:?}/fcsgx_stage_\${sample}"
    mkdir -p "\$STAGE"
    rsync -a --delete ${params.fcsgx_db}/ "\$STAGE/"
    if [ ! -s "\$STAGE/all.gxi" ]; then
        echo "[ERROR] FCS-GX DB staging failed: \$STAGE/all.gxi missing (is ${params.fcsgx_db} bound into the container?)" >&2
        exit 1
    fi
    AAFTF fcs_gx_purge --db "\$STAGE/all" \\
        -t "${phylum}" \\
        -i ${asm} -o ${sample}.contam_clean.fasta \\
        -w fcsgx_work_${sample}
    rm -rf "\$STAGE"
    gzip -c ${sample}.contam_clean.fasta > ${sample}.contam_clean.fasta.gz
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'fcs_gx=command -v gx >/dev/null && ls \$(dirname \$(command -v gx))/../conda-meta 2>/dev/null | grep ncbi-fcs-gx | head -1' \\
        "fcs_gx_db=readlink -f ${params.fcsgx_db}"
    """

    stub:
    """
    cp ${asm} ${sample}.contam_clean.fasta
    gzip -c ${sample}.contam_clean.fasta > ${sample}.contam_clean.fasta.gz
    """
}

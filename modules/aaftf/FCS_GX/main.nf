process FCS_GX {
    tag   { sample }
    label 'aaftf'
    publishDir "${params.outdir}/fcs_gx", mode: 'copy', pattern: '*.fcs_gx.fasta'

    input:
    tuple val(sample), path(asm), val(taxonid)

    output:
    tuple val(sample), path("${sample}.fcs_gx.fasta"), emit: clean

    script:
    """
    # Optional NCBI FCS-GX contamination purge. taxonid is the NCBI taxonomy id
    # of the ORGANISM (e.g. 4890 for Ascomycota / 4751 for Fungi) taken from the
    # samples.csv `taxid` column (or params.fcs_taxid fallback).
    #
    # FCS-GX -t expects a taxid at the PHYLUM level. Resolve it at runtime from
    # the organism id with taxonkit so hardcoding drift is avoided (mirrors the
    # reference BFD GENOME_CLEAN approach): lineage -> "{p}" phylum name ->
    # name2taxid -> numeric phylum taxid.
    TAXONKIT_DB=${params.taxondb}
    phylum=\$(echo "${taxonid}" | taxonkit --data-dir \$TAXONKIT_DB lineage \\
        | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result \\
        | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
    if [ -z "\$phylum" ]; then
        echo "FCS_GX: could not resolve phylum taxid from ${taxonid}, falling back to params.fcs_taxid (${params.fcs_taxid})" >&2
        phylum=${params.fcs_taxid}
    fi
    if ! [[ "\$phylum" =~ ^[0-9]+\$ ]]; then
        echo "FCS_GX: ERROR: resolved phylum '\$phylum' is not numeric; cannot run fcs_gx_purge for ${sample}" >&2
        exit 1
    fi
    echo "FCS_GX: ${sample} taxonid=${taxonid} -> phylum_taxid=\$phylum"

    # The FCS-GX database is large; stage it to node-local scratch once per sample.
    STAGE="\${SCRATCH:?}/fcsgx_stage_\${sample}"
    mkdir -p "\$STAGE"
    rsync -a --delete ${params.fcsgx_db}/ "\$STAGE/"
    AAFTF fcs_gx_purge --db "\$STAGE/all" \\
        -t "\$phylum" -c ${task.cpus} \\
        -i ${asm} -o ${sample}.fcs_gx.fasta \\
        -w fcsgx_work_${sample}
    rm -rf "\$STAGE"
    """

    stub:
    """
    cp ${asm} ${sample}.fcs_gx.fasta
    """
}

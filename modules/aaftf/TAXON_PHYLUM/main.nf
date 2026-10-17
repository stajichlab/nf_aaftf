// Resolve the phylum of a sample's organism taxid with taxonkit.
// taxonkit runs in its own container (params.taxonkit_container) because the
// AAFTF >= 0.7 image no longer ships it.
//   phylum_taxid.txt  numeric phylum taxid for FCS-GX (-t). Falls back to
//                     params.fcs_taxid when the lookup returns nothing.
//   phylum_name.txt   phylum name for `AAFTF taxpurge -p`. Empty if unresolved.
process TAXON_PHYLUM {
    tag   { sample }
    label 'taxonkit'

    input:
    tuple val(sample), val(taxonid)

    output:
    tuple val(sample), path('phylum_taxid.txt'), path('phylum_name.txt'), emit: phylum

    script:
    """
    TAXONKIT_DB=${params.taxondb}
    name=\$(echo "${taxonid}" | taxonkit --data-dir \$TAXONKIT_DB lineage \\
        | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result \\
        | cut -f3 | head -n 1 || true)
    if [ "\$name" = "No rank" ]; then name=""; fi
    taxid=\$(echo "${taxonid}" | taxonkit --data-dir \$TAXONKIT_DB lineage \\
        | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result \\
        | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1 || true)
    if [ -z "\$taxid" ]; then
        echo "TAXON_PHYLUM: could not resolve phylum taxid from ${taxonid} for ${sample}; using params.fcs_taxid (${params.fcs_taxid})" >&2
        taxid=${params.fcs_taxid}
    fi
    echo "TAXON_PHYLUM: ${sample} taxonid=${taxonid} -> phylum='\$name' phylum_taxid=\$taxid"
    printf '%s\\n' "\$taxid" > phylum_taxid.txt
    printf '%s\\n' "\$name"  > phylum_name.txt
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'taxonkit=taxonkit version'
    """

    stub:
    """
    echo ${params.fcs_taxid} > phylum_taxid.txt
    echo Ascomycota > phylum_name.txt
    """
}

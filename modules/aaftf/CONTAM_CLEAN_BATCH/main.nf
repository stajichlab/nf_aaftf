// Batched FCS-GX purge for a list of assemblies.
//
// The FCS-GX database is huge (~465 GB), so rsync-staging it once per SAMPLE
// (what CONTAM_CLEAN does) pays that full sync cost on every genome. This
// process stages it ONCE per batch into node-local $SCRATCH and then loops
// AAFTF fcs_gx_purge across every assembly in the batch, amortizing the sync
// over up to contam_clean_batch_size genomes (each purge itself only takes
// ~1-2 min once the DB is staged).
//
// Cleaned assemblies are written directly to a fixed location
// (params.outdir/contam_clean/<sample>.contam_clean.fasta, matching what
// CONTAM_CLEAN's own publishDir produces) instead of through Nextflow's
// per-task output staging, and skipped here if already present. That means
// a batch that fails partway, or a whole pipeline re-launch, does not repay
// the DB-staging cost for genomes a prior attempt already cleaned — main.nf
// filters those out before batching and re-derives the per-sample channel
// from this fixed path afterwards. Mirrors GENOME_CLEAN_BATCH /
// FUNANNOTATE_GENOME_PREP in nf_funannotate1.
process CONTAM_CLEAN_BATCH {
    tag { "batch_${task.index}" }
    label 'aaftf'

    input:
    tuple val(items), path(assemblies)

    output:
    path "contam_clean_batch_${task.index}.manifest.tsv", emit: manifest

    script:
    // items: List of [sample, filename, taxonid], built in main.nf. `filename`
    // is the basename Nextflow stages the matching entry of `assemblies`
    // under in this task's work dir, so there is no ordering ambiguity
    // between the metadata list and the staged files.
    def batch_tsv = items.collect { s, fname, taxonid -> "${s}\t${fname}\t${taxonid}" }.join('\n')
    """
    set -uo pipefail
    TAXONKIT_DB=${params.taxondb}
    DEST=${params.outdir}/contam_clean
    mkdir -p \$DEST

    cat > batch.tsv <<'BATCH_EOF'
${batch_tsv}
BATCH_EOF

    MANIFEST=contam_clean_batch_${task.index}.manifest.tsv
    : > \$MANIFEST
    n_total=\$(grep -c . batch.tsv || true)
    echo "[INFO] batch ${task.index}: \$n_total genomes to consider"

    # Stage the FCS-GX DB into node-local scratch ONCE for the whole batch.
    STAGE="\${SCRATCH:?}/fcsgx_stage_batch_${task.index}"
    mkdir -p "\$STAGE"
    rsync -a --delete ${params.fcsgx_db}/ "\$STAGE/"

    i=0
    fcs_fails=0
    while IFS=\$'\\t' read -r sample fname taxonid; do
        [ -z "\$sample" ] && continue
        i=\$((i+1))
        target=\$DEST/\${sample}.contam_clean.fasta
        if [ -s "\$target" ]; then
            echo "[\$i/\$n_total][SKIP] \$sample already cleaned"
            printf '%s\\t%s\\n' "\$sample" "\$target" >> \$MANIFEST
            continue
        fi

        phylum=\$(echo "\$taxonid" | taxonkit --data-dir \$TAXONKIT_DB lineage \\
            | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result \\
            | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
        if [ -z "\$phylum" ]; then
            echo "[\$i/\$n_total][WARN] could not resolve phylum for \$sample (taxonid=\$taxonid), falling back to params.fcs_taxid (${params.fcs_taxid})" >&2
            phylum=${params.fcs_taxid}
        fi
        if ! [[ "\$phylum" =~ ^[0-9]+\$ ]]; then
            echo "[\$i/\$n_total][FAIL] resolved phylum '\$phylum' is not numeric for \$sample" >&2
            fcs_fails=\$((fcs_fails+1))
            continue
        fi
        echo "[\$i/\$n_total][INFO] \$sample taxonid=\$taxonid -> phylum_taxid=\$phylum"

        if AAFTF fcs_gx_purge --db "\$STAGE/all" \\
            -t "\$phylum" -c ${task.cpus} \\
            -i "\$fname" -o "\${target}.tmp" \\
            -w fcsgx_work_\${sample}; then
            mv "\${target}.tmp" "\$target"
            echo "[\$i/\$n_total][OK] \$sample -> \$target"
            printf '%s\\t%s\\n' "\$sample" "\$target" >> \$MANIFEST
        else
            echo "[\$i/\$n_total][FAIL] fcs_gx_purge failed for \$sample" >&2
            fcs_fails=\$((fcs_fails+1))
        fi
    done < batch.tsv
    rm -rf "\$STAGE"

    n_done=\$(grep -c . \$MANIFEST || echo 0)
    if [ "\$fcs_fails" -gt 0 ] || [ "\$n_done" -eq 0 ]; then
        echo "[ERROR] batch ${task.index}: \$fcs_fails purge failure(s); manifest has \$n_done entries" >&2
        exit 1
    fi
    echo "[INFO] batch ${task.index} complete: \$n_done cleaned assemblies"
    """

    stub:
    def batch_tsv = items.collect { s, fname, taxonid -> "${s}\t${fname}\t${taxonid}" }.join('\n')
    """
    DEST=${params.outdir}/contam_clean
    mkdir -p \$DEST
    MANIFEST=contam_clean_batch_${task.index}.manifest.tsv
    : > \$MANIFEST
    cat > batch.tsv <<'BATCH_EOF'
${batch_tsv}
BATCH_EOF
    while IFS=\$'\\t' read -r sample fname taxonid; do
        [ -z "\$sample" ] && continue
        cp "\$fname" \$DEST/\${sample}.contam_clean.fasta
        printf '%s\\t%s\\n' "\$sample" "\$DEST/\${sample}.contam_clean.fasta" >> \$MANIFEST
    done < batch.tsv
    """
}

// Batched FCS-GX purge for a list of assemblies.
//
// The FCS-GX database is huge (~465 GB), so rsync-staging it once per SAMPLE
// (what CONTAM_CLEAN does) pays that full sync cost on every genome. This
// process stages it ONCE per batch into ${params.fcsgx_stage_dir} (default /dev/shm) and then loops
// AAFTF fcs_gx_purge across every assembly in the batch, amortizing the sync
// over up to contam_clean_batch_size genomes (each purge itself only takes
// ~1-2 min once the DB is staged).
//
// Cleaned assemblies are written directly to a fixed location
// (params.outdir/contam_clean/<sample>.contam_clean.fasta.gz, matching what
// CONTAM_CLEAN's own publishDir produces, gzip-compressed to save space)
// instead of through Nextflow's per-task output staging, and skipped here if
// already present. That means a batch that fails partway, or a whole
// pipeline re-launch, does not repay the DB-staging cost for genomes a prior
// attempt already cleaned — main.nf filters those out before batching and
// re-derives the per-sample channel from this fixed path afterwards (via a
// GUNZIP step, since fcs_gx_purge itself needs a plain FASTA). Mirrors
// GENOME_CLEAN_BATCH / FUNANNOTATE_GENOME_PREP in nf_funannotate1.
process CONTAM_CLEAN_BATCH {
    tag { "batch_${task.index}" }
    label 'aaftf'

    input:
    tuple val(items), path(assemblies)

    output:
    path "contam_clean_batch_${task.index}.manifest.tsv", emit: manifest

    script:
    // items: List of [sample, filename, phylum_taxid], built in main.nf. `filename`
    // is the basename Nextflow stages the matching entry of `assemblies`
    // under in this task's work dir, so there is no ordering ambiguity
    // between the metadata list and the staged files.
    def batch_tsv = items.collect { s, fname, phylum -> "${s}\t${fname}\t${phylum}" }.join('\n')
    """
    set -uo pipefail
    DEST=${params.outdir}/contam_clean
    mkdir -p \$DEST

    cat > batch.tsv <<'BATCH_EOF'
${batch_tsv}
BATCH_EOF

    MANIFEST=contam_clean_batch_${task.index}.manifest.tsv
    : > \$MANIFEST
    n_total=\$(grep -c . batch.tsv || true)
    echo "[INFO] batch ${task.index}: \$n_total genomes to consider"

    # Stage the FCS-GX DB into ${params.fcsgx_stage_dir} (default /dev/shm, RAM-backed) ONCE for the whole batch.
    # The trap removes the staged copy on any exit, including failure and scancel (SIGTERM).
    # Nextflow starts the container with `env -`, so USER and SLURM_JOB_ID are unset (and set -u would stop on them).
    # mktemp -d gives each task its own directory, so two tasks on one node cannot share or delete each other's copy.
    STAGE="\$(mktemp -d "${params.fcsgx_stage_dir}/fcsgx_batch_${task.index}_XXXXXX")" || { echo "[ERROR] cannot create a stage dir under ${params.fcsgx_stage_dir}" >&2; exit 1; }
    trap 'rm -rf "\$STAGE"' EXIT
    trap 'exit 143' TERM INT
    mkdir -p "\$STAGE"
    need_kb=\$(du -sk --apparent-size ${params.fcsgx_db} | cut -f1)
    have_kb=\$(df -k --output=avail "\$STAGE" | tail -1)
    if [ "\$have_kb" -lt "\$need_kb" ]; then
        echo "[ERROR] \$STAGE has \$have_kb KB free; FCS-GX DB needs \$need_kb KB. Set --fcsgx_stage_dir to a larger filesystem." >&2
        exit 1
    fi
    echo "[INFO] stage dir \$STAGE: \$(df -h --output=size,avail "\$STAGE" | tail -1 | xargs) (size avail); cgroup memory.max: \$(cat /sys/fs/cgroup/memory.max 2>/dev/null || echo n/a)"
    echo "[TIME] stage_start \$(date +%s)"
    rsync -a --delete ${params.fcsgx_db}/ "\$STAGE/"
    echo "[TIME] stage_end \$(date +%s)"
    if [ ! -s "\$STAGE/all.gxi" ]; then
        echo "[ERROR] FCS-GX DB staging failed: \$STAGE/all.gxi missing (is ${params.fcsgx_db} bound into the container?)" >&2
        exit 1
    fi

    i=0
    fcs_fails=0
    while IFS=\$'\\t' read -r sample fname phylum; do
        [ -z "\$sample" ] && continue
        i=\$((i+1))
        target=\$DEST/\${sample}.contam_clean.fasta.gz
        if [ -s "\$target" ]; then
            echo "[\$i/\$n_total][SKIP] \$sample already cleaned"
            printf '%s\\t%s\\n' "\$sample" "\$target" >> \$MANIFEST
            continue
        fi

        if ! [[ "\$phylum" =~ ^[0-9]+\$ ]]; then
            echo "[\$i/\$n_total][FAIL] resolved phylum '\$phylum' is not numeric for \$sample" >&2
            fcs_fails=\$((fcs_fails+1))
            continue
        fi
        echo "[\$i/\$n_total][INFO] \$sample phylum_taxid=\$phylum"

        # AAFTF >= 0.7: -o is an output DIRECTORY that receives <input prefix>.fcsgx_cleaned.fasta and
        # <input prefix>.fcs_gx_report.txt. Link the input as <sample>.fa so the prefix is the sample.
        mkdir -p fcs_in "fcs_out_\${sample}"
        ln -sf "\$(readlink -f "\$fname")" "fcs_in/\${sample}.fa"
        cleaned_tmp="fcs_out_\${sample}/\${sample}.fcsgx_cleaned.fasta"
        # The exit code is not a reliable success signal; require the cleaned file.
        if AAFTF fcs_gx_purge -d "\$STAGE/all" \\
            -t "\$phylum" \\
            -i "fcs_in/\${sample}.fa" -c ${task.cpus} -o "fcs_out_\${sample}" \\
            -w fcsgx_work_\${sample} \\
            && [ -s "\$cleaned_tmp" ]; then
            gzip -c "\$cleaned_tmp" > "\${target}.tmp"
            mv "\${target}.tmp" "\$target"
            rm -rf "fcs_out_\${sample}"
            echo "[\$i/\$n_total][OK] \$sample -> \$target"
            printf '%s\\t%s\\n' "\$sample" "\$target" >> \$MANIFEST
        else
            echo "[\$i/\$n_total][FAIL] fcs_gx_purge failed for \$sample" >&2
            rm -rf "fcs_out_\${sample}"
            fcs_fails=\$((fcs_fails+1))
        fi
    done < batch.tsv

    n_done=\$(grep -c . \$MANIFEST || echo 0)
    if [ "\$fcs_fails" -gt 0 ] || [ "\$n_done" -eq 0 ]; then
        echo "[ERROR] batch ${task.index}: \$fcs_fails purge failure(s); manifest has \$n_done entries" >&2
        exit 1
    fi
    echo "[INFO] batch ${task.index} complete: \$n_done cleaned assemblies"
    while IFS=\$'\\t' read -r s_id _ _; do
        [ -z "\$s_id" ] && continue
        ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/\${s_id}/versions" "${task.process}" "\$s_id" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
            'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
            'fcs_gx=command -v gx >/dev/null && ls \$(dirname \$(command -v gx))/../conda-meta 2>/dev/null | grep ncbi-fcs-gx | head -1' \\
            "fcs_gx_db=readlink -f ${params.fcsgx_db}"
    done < batch.tsv
    """

    stub:
    def batch_tsv = items.collect { s, fname, phylum -> "${s}\t${fname}\t${phylum}" }.join('\n')
    """
    DEST=${params.outdir}/contam_clean
    mkdir -p \$DEST
    MANIFEST=contam_clean_batch_${task.index}.manifest.tsv
    : > \$MANIFEST
    cat > batch.tsv <<'BATCH_EOF'
${batch_tsv}
BATCH_EOF
    while IFS=\$'\\t' read -r sample fname phylum; do
        [ -z "\$sample" ] && continue
        gzip -c "\$fname" > \$DEST/\${sample}.contam_clean.fasta.gz
        printf '%s\\t%s\\n' "\$sample" "\$DEST/\${sample}.contam_clean.fasta.gz" >> \$MANIFEST
    done < batch.tsv
    """
}

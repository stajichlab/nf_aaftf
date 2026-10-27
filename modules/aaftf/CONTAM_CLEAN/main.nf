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

    # The FCS-GX database is large; stage it into ${params.fcsgx_stage_dir} (default /dev/shm, RAM-backed) once per sample.
    # The trap removes the staged copy on any exit, including failure and scancel (SIGTERM).
    STAGE="${params.fcsgx_stage_dir}/fcsgx_\${USER}_\${SLURM_JOB_ID:-\$\$}_${sample}"
    trap 'rm -rf "\$STAGE"' EXIT
    trap 'exit 143' TERM INT
    mkdir -p "\$STAGE"
    need_kb=\$(du -sk --apparent-size ${params.fcsgx_db} | cut -f1)
    have_kb=\$(df -k --output=avail "\$STAGE" | tail -1)
    if [ "\$have_kb" -lt "\$need_kb" ]; then
        echo "[ERROR] \$STAGE has \$have_kb KB free; FCS-GX DB needs \$need_kb KB. Set --fcsgx_stage_dir to a larger filesystem." >&2
        exit 1
    fi
    rsync -a --delete ${params.fcsgx_db}/ "\$STAGE/"
    if [ ! -s "\$STAGE/all.gxi" ]; then
        echo "[ERROR] FCS-GX DB staging failed: \$STAGE/all.gxi missing (is ${params.fcsgx_db} bound into the container?)" >&2
        exit 1
    fi
    # AAFTF >= 0.7: -o is an output DIRECTORY that receives <input prefix>.fcsgx_cleaned.fasta and
    # <input prefix>.fcs_gx_report.txt. Link the input as <sample>.fa so the prefix is the sample.
    mkdir -p fcs_in fcs_out
    ln -sf "\$(readlink -f ${asm})" fcs_in/${sample}.fa
    AAFTF fcs_gx_purge -d "\$STAGE/all" \\
        -t "${phylum}" \\
        -i fcs_in/${sample}.fa -c ${task.cpus} -o fcs_out \\
        -w fcsgx_work_${sample}
    # The exit code is not a reliable success signal; require the cleaned file.
    if [ ! -s fcs_out/${sample}.fcsgx_cleaned.fasta ]; then
        echo "CONTAM_CLEAN: ERROR: fcs_gx_purge produced no cleaned FASTA for ${sample}" >&2
        exit 1
    fi
    cp fcs_out/${sample}.fcsgx_cleaned.fasta ${sample}.contam_clean.fasta
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

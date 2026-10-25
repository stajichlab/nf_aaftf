process FILTER {
    tag   { sample }
    label 'aaftf_trim'
    // Previous pattern was not a valid Nextflow glob (see TRIM/main.nf), so
    // nothing was ever actually copied to outdir/filter -- confirmed empty
    // 2026-09-26. No pattern = copy exactly the declared `output` paths.
    publishDir "${params.outdir}/filter", mode: 'copy'

    input:
    tuple val(sample), path(reads1), path(reads2), path(merged)

    output:
    tuple val(sample), path("${sample}_filtered_1.fastq.gz"), path("${sample}_filtered_2.fastq.gz"), path("${sample}_filtered_U.fastq.gz"), emit: filtered

    script:
    """
    # AAFTF filter does not check the exit status of shuffle.sh, bbduk.sh or reformat.sh.
    # If one is killed (SLURM OUT_OF_MEMORY), the output is missing and AAFTF fails later
    # ("pigz -dc ..._filtered_1.fastq.gz returned non-zero exit status 1", AAFTF exit 1).
    # Without -e the script ignored that exit status, ran the single-end filter and exited 0.
    # Stop at the first failure and check each output so Nextflow retries with more memory.
    set -euo pipefail
    trap 'echo "FILTER: ERROR: command failed at line \$LINENO of .command.sh (exit \$?)" >&2' ERR
    export AAFTF_DB=/opt/aaftf_db   # AAFTF >= 0.7 reads the database folder from the AAFTF_DB environment variable (the --AAFTF_DB flag was removed)
    # Paired-end reads: filter against phiX + common contaminant/vector DBs.
    AAFTF filter --aligner ${params.filter_aligner} \
        -c ${task.cpus} -m 16 \
        -1 ${reads1} -2 ${reads2} \
        -o ${sample}
    test -s ${sample}_filtered_1.fastq.gz && test -s ${sample}_filtered_2.fastq.gz

    # Merged (single-end) reads: same filter run as single-end so merged reads
    # survive too. AAFTF >= 0.7 names single-end output <prefix>_filtered_1 (it used to
    # be _filtered_U), which would overwrite the paired _filtered_1 above. Use a
    # separate prefix, then rename to the _filtered_U the pipeline expects.
    AAFTF filter --aligner ${params.filter_aligner} \
        -c ${task.cpus} -m 16 \
        -1 ${merged} \
        -o ${sample}_SE
    test -s ${sample}_SE_filtered_1.fastq.gz
    mv ${sample}_SE_filtered_1.fastq.gz ${sample}_filtered_U.fastq.gz
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'bbduk=bbduk.sh --version 2>&1 | grep -i '\\''BBTools version'\\''' \\
        'bwa=bwa 2>&1 | grep -i '\\''^Version'\\''' \\
        'minimap2=minimap2 --version' \\
        'samtools=samtools --version | head -1' \\
        "filter_aligner=echo ${params.filter_aligner}"
    """

    stub:
    """
    touch ${sample}_filtered_1.fastq.gz
    touch ${sample}_filtered_2.fastq.gz
    touch ${sample}_filtered_U.fastq.gz
    """
}

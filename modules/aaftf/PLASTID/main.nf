process PLASTID {
    tag   { sample }
    label 'aaftf_lite'
    publishDir "${params.outdir}/plastid", mode: 'copy',
        pattern: '*.plastid.{fasta.gz,validation.tsv}'

    input:
    tuple val(sample), path(asm), path(depth_report)

    output:
    tuple val(sample), path("${sample}.noplastid.fasta"), emit: nuclear
    path "${sample}.plastid.fasta.gz"
    path "${sample}.plastid.validation.tsv"

    script:
    def veto_cov = params.getOrDefault('rdna_veto_cov', 0.2)
    """
    set -euo pipefail
    trap 'echo "PLASTID: ERROR: command failed at line \$LINENO of .command.sh (exit \$?)" >&2' ERR
    # (the aaftf profile sets shell=/bin/bash without -ue, so fail fast explicitly)
    DB=${params.organelle_db}
    VALIDATE=${projectDir}/bin/organelle_validate.py
    # Chloroplast/plastid contigs are high copy but not AT-rich, so AAFTF depth's MITO rule does
    # not see them. Candidates are the depth OUTLIER and ELEVATED contigs of the MITOCHONDRION depth
    # report that are still in the assembly (mito contigs are already gone). Only candidates are
    # searched against RefSeq plastid genomes; a contig is removed only when alignments cover at
    # least plastid_min_cov of its length and it is not rDNA (>= rdna_veto_cov of fungal 18S/28S/ITS).
    if [ ! -s "\$DB/plastid/plastid.nhr" ]; then
        echo "PLASTID: ERROR: no plastid BLAST DB at \$DB/plastid; run bin/download_organelle_db.sh plastid" >&2
        exit 1
    fi
    awk '/\\*\\* OUTLIER \\(|ELEVATED \\(inspect/ {print \$1}' ${depth_report} > flagged.ids
    awk -v ids=flagged.ids 'BEGIN{while((getline l < ids)>0) k[l]=1} /^>/{id=substr(\$1,2); p=(id in k)} p' ${asm} > candidates.fasta
    \$VALIDATE --candidates candidates.fasta --db \$DB/plastid/plastid \\
        --min-cov ${params.plastid_min_cov} --min-pident ${params.organelle_min_pident} --cpus ${task.cpus} \\
        --out plastid_val.tsv --validated plastid_val.ids

    : > rdna.tsv
    if [ -s candidates.fasta ] && [ "${veto_cov}" != "0" ]; then
        if [ ! -s "\$DB/rdna/rdna.nal" ]; then
            echo "PLASTID: ERROR: no rDNA BLAST DB at \$DB/rdna; run bin/download_organelle_db.sh rdna (or set --rdna_veto_cov 0)" >&2
            exit 1
        fi
        \$VALIDATE --candidates candidates.fasta --db \$DB/rdna/rdna --min-cov 2 --min-pident 80 \\
            --cpus ${task.cpus} --out rdna.tsv --validated /dev/null
    fi
    awk -F'\\t' -v OFS='\\t' -v veto=${veto_cov} '
        FILENAME=="rdna.tsv" { if (FNR>1) rd[\$1]=\$4; next }
        FNR==1 { print \$0, "rdna_frac", "decision"; next }
        { r = (\$1 in rd) ? rd[\$1] : 0
          d = (\$8!="yes") ? "kept_unvalidated" : ((veto>0 && r>=veto) ? "kept_rdna" : "removed")
          print \$0, r, d }' rdna.tsv plastid_val.tsv > ${sample}.plastid.validation.tsv
    awk -F'\\t' 'NR>1 && \$NF=="removed" {print \$1}' ${sample}.plastid.validation.tsv > remove.ids

    awk -v ids=remove.ids 'BEGIN{while((getline l < ids)>0) k[l]=1} /^>/{id=substr(\$1,2); p=(id in k)} p' ${asm} > ${sample}.plastid.fasta
    awk -v ids=remove.ids 'BEGIN{while((getline l < ids)>0) k[l]=1} /^>/{id=substr(\$1,2); p=!(id in k)} p' ${asm} > ${sample}.noplastid.fasta
    n_in=\$(grep -c '^>' ${asm} || true); n_keep=\$(grep -c '^>' ${sample}.noplastid.fasta || true); n_rm=\$(grep -c '^>' ${sample}.plastid.fasta || true)
    if [ "\$((n_keep + n_rm))" -ne "\$n_in" ] || [ "\$n_keep" -eq 0 ]; then
        echo "PLASTID: ERROR: split mismatch for ${sample}: in=\$n_in kept=\$n_keep plastid=\$n_rm" >&2
        exit 1
    fi
    echo "PLASTID: ${sample}: \$n_in contigs -> \$n_keep kept, \$n_rm plastid removed (\$(grep -c . flagged.ids || true) candidates)"
    gzip -c ${sample}.plastid.fasta > ${sample}.plastid.fasta.gz
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'blastn=blastn -version | head -1' \\
        "plastid_db=cat ${params.organelle_db}/plastid/RELEASE_NUMBER 2>/dev/null || echo unknown"
    """

    stub:
    """
    cp ${asm} ${sample}.noplastid.fasta
    : > ${sample}.plastid.fasta
    gzip -c ${sample}.plastid.fasta > ${sample}.plastid.fasta.gz
    printf 'contig\\tdecision\\n' > ${sample}.plastid.validation.tsv
    """
}

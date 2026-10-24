process MITOCHONDRION {
    tag   { sample }
    label 'aaftf_lite'
    publishDir "${params.outdir}/mitochondrion", mode: 'copy',
        pattern: '*.mitochondrion.{fasta.gz,candidates.tsv,validation.tsv,depth.txt,warnings.txt}'

    input:
    tuple val(sample), path(asm), path(filtered_1), path(filtered_2)

    output:
    tuple val(sample), path("${sample}.nuclear.fasta"), path("${sample}.mitochondrion.depth.txt"), emit: nuclear
    tuple val(sample), path("${sample}.mitochondrion.warnings.txt"), emit: warnings
    path "${sample}.mitochondrion.fasta.gz"
    path "${sample}.mitochondrion.candidates.tsv"
    path "${sample}.mitochondrion.validation.tsv"

    script:
    def require_hit = (params.getOrDefault('mitochondrion_require_hit', true) as String).toBoolean()
    def veto_cov    = params.getOrDefault('rdna_veto_cov', 0.2)
    """
    set -euo pipefail
    trap 'echo "MITOCHONDRION: ERROR: command failed at line \$LINENO of .command.sh (exit \$?)" >&2' ERR
    # (the aaftf profile sets shell=/bin/bash without -ue, so fail fast explicitly)
    # Fail early with a clear message if this AAFTF image predates the depth MITO flag.
    depth_help=\$(AAFTF depth -h 2>&1 || true)
    if ! grep -q -- '--mito_fold' <<< "\$depth_help"; then
        echo "MITOCHONDRION: ERROR: this AAFTF image has no 'depth --mito_fold' option; use an image built from AAFTF main >= 0904010 (set --aaftf_sif)" >&2
        exit 1
    fi
    DB=${params.organelle_db}
    VALIDATE=${projectDir}/bin/organelle_validate.py
    # fasta_of IDS FASTA > subset
    fasta_of() { awk -v ids="\$1" 'BEGIN{while((getline l < ids)>0) k[l]=1} /^>/{id=substr(\$1,2); p=(id in k)} p' "\$2"; }

    # 1. Map the filtered reads back to the contamination-cleaned assembly (before RMDUP) and flag
    #    candidate mitochondrial contigs (AAFTF depth: depth >= mito_fold x nuclear baseline,
    #    GC >= mito_gc_drop points below the assembly, length >= mito_min_len). Needs an AAFTF with
    #    the MITO flag (commit 0904010 or later; not in the 0.7.0-beta.4-18 SIF).
    AAFTF depth -i ${asm} -o ${sample}.mitochondrion.depth.txt \\
        -1 ${filtered_1} -2 ${filtered_2} --no-plot \\
        --mito_fold ${params.mito_fold} --mito_gc_drop ${params.mito_gc_drop} --mito_min_len ${params.mito_min_len} \\
        -c ${task.cpus} -w depth_work_${sample}

    TSV=${sample}.mitochondrion.depth.mito_candidates.tsv
    if [ ! -s "\$TSV" ]; then
        echo "MITOCHONDRION: ERROR: AAFTF depth wrote no \$TSV; the container needs the MITO flag (AAFTF main >= 0904010)" >&2
        exit 1
    fi
    mv "\$TSV" ${sample}.mitochondrion.candidates.tsv
    awk -F'\\t' 'NR>1 && \$6=="yes" {print \$1}' ${sample}.mitochondrion.candidates.tsv > candidates.ids
    fasta_of candidates.ids ${asm} > candidates.fasta

    # 2. Validate ONLY the candidates against RefSeq mitochondrial genomes by whole-contig coverage
    #    (union of alignments / contig length), so a NUMT insertion in a nuclear contig does not count.
    if [ "${require_hit}" = "true" ]; then
        if [ ! -s "\$DB/mitochondrion/mitochondrion.nhr" ]; then
            echo "MITOCHONDRION: ERROR: no mitochondrion BLAST DB at \$DB/mitochondrion; run bin/download_organelle_db.sh (or set --mitochondrion_require_hit false)" >&2
            exit 1
        fi
        \$VALIDATE --candidates candidates.fasta --db \$DB/mitochondrion/mitochondrion \\
            --min-cov ${params.mito_min_cov} --min-pident ${params.organelle_min_pident} --cpus ${task.cpus} \\
            --out mito_val.tsv --validated mito_val.ids
    else
        printf 'contig\\tlength\\tcovered_bp\\tcovered_frac\\tbest_pident\\tn_hsps\\ttop_subject\\tvalidated\\n' > mito_val.tsv
        awk '{print \$1"\\t-\\t-\\t-\\t-\\t-\\t-\\tyes"}' candidates.ids >> mito_val.tsv
    fi

    # 3. rDNA veto: nuclear rDNA is high copy and can look like an organelle, so a contig that
    #    is >= rdna_veto_cov (${veto_cov}) rDNA (NCBI fungal 18S/28S/ITS) is never removed.
    #    The same screen labels depth OUTLIER contigs in the warnings (rDNA-like or unexplained).
    awk '/\\*\\* OUTLIER/ && !/MITO/ {print \$1}' ${sample}.mitochondrion.depth.txt > outlier.ids
    cat candidates.ids outlier.ids | sort -u > rdna_query.ids
    : > rdna.tsv
    if [ -s rdna_query.ids ] && [ "${veto_cov}" != "0" ]; then
        if [ ! -s "\$DB/rdna/rdna.nal" ]; then
            echo "MITOCHONDRION: ERROR: no rDNA BLAST DB at \$DB/rdna; run bin/download_organelle_db.sh rdna (or set --rdna_veto_cov 0)" >&2
            exit 1
        fi
        fasta_of rdna_query.ids ${asm} > rdna_query.fasta
        \$VALIDATE --candidates rdna_query.fasta --db \$DB/rdna/rdna --min-cov 2 --min-pident 80 \\
            --cpus ${task.cpus} --out rdna.tsv --validated /dev/null
    fi

    # Final table for the candidates: validator columns + rDNA coverage + decision.
    awk -F'\\t' -v OFS='\\t' -v veto=${veto_cov} '
        FILENAME=="rdna.tsv" { if (FNR>1) rd[\$1]=\$4; next }
        FNR==1 { print \$0, "rdna_frac", "decision"; next }
        { r = (\$1 in rd) ? rd[\$1] : 0
          d = (\$8!="yes") ? "kept_unvalidated" : ((veto>0 && r>=veto) ? "kept_rdna" : "removed")
          print \$0, r, d }' rdna.tsv mito_val.tsv > ${sample}.mitochondrion.validation.tsv
    awk -F'\\t' 'NR>1 && \$NF=="removed" {print \$1}' ${sample}.mitochondrion.validation.tsv > remove.ids

    # 4. Split: contigs in remove.ids go to the mitochondrion file, all others continue.
    awk -v ids=remove.ids 'BEGIN{while((getline l < ids)>0) k[l]=1} /^>/{id=substr(\$1,2); p=(id in k)} p' ${asm} > ${sample}.mitochondrion.fasta
    awk -v ids=remove.ids 'BEGIN{while((getline l < ids)>0) k[l]=1} /^>/{id=substr(\$1,2); p=!(id in k)} p' ${asm} > ${sample}.nuclear.fasta
    n_in=\$(grep -c '^>' ${asm} || true); n_nuc=\$(grep -c '^>' ${sample}.nuclear.fasta || true); n_mt=\$(grep -c '^>' ${sample}.mitochondrion.fasta || true)
    if [ "\$((n_nuc + n_mt))" -ne "\$n_in" ] || [ "\$n_nuc" -eq 0 ]; then
        echo "MITOCHONDRION: ERROR: split mismatch for ${sample}: in=\$n_in nuclear=\$n_nuc mitochondrion=\$n_mt" >&2
        exit 1
    fi
    echo "MITOCHONDRION: ${sample}: \$n_in contigs -> \$n_nuc kept, \$n_mt mitochondrial removed (\$(wc -l < candidates.ids) candidates)"
    gzip -c ${sample}.mitochondrion.fasta > ${sample}.mitochondrion.fasta.gz

    # 5. Warnings (nothing removed for these): MITO candidates that were kept, and depth OUTLIER
    #    contigs that are not MITO, labeled by rDNA coverage.
    awk -F'\\t' 'NR>1 && \$NF!="removed" {print "MITO_" toupper(\$NF) "\\t" \$1 "\\tlen=" \$2 "\\tmito_cov=" \$4 "\\trdna_cov=" \$(NF-1)}' ${sample}.mitochondrion.validation.tsv > ${sample}.mitochondrion.warnings.txt
    awk -v OFS='\\t' -v veto=${veto_cov} '
        FILENAME=="rdna.tsv" { if (FNR>1) rd[\$1]=\$4; next }
        /\\*\\* OUTLIER/ && !/MITO/ { r = (\$1 in rd) ? rd[\$1] : 0
            print "OUTLIER", \$1, "len=" \$2, "depth=" \$3, "GC=" \$4, "rdna_cov=" r, (r>=veto ? "rDNA-like" : "unexplained") }' rdna.tsv ${sample}.mitochondrion.depth.txt >> ${sample}.mitochondrion.warnings.txt
    if [ -s ${sample}.mitochondrion.warnings.txt ]; then
        echo "MITOCHONDRION: WARNING: ${sample}: contigs kept that need review:" >&2
        cat ${sample}.mitochondrion.warnings.txt >&2
    fi
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'minimap2=minimap2 --version' \\
        'samtools=samtools --version | head -1' \\
        'mosdepth=mosdepth --version' \\
        'blastn=blastn -version | head -1' \\
        "mitochondrion_db=cat ${params.organelle_db}/mitochondrion/RELEASE_NUMBER 2>/dev/null || echo unknown"
    """

    stub:
    """
    cp ${asm} ${sample}.nuclear.fasta
    : > ${sample}.mitochondrion.fasta
    gzip -c ${sample}.mitochondrion.fasta > ${sample}.mitochondrion.fasta.gz
    printf 'contig\\tlength\\tmean_depth\\tfold_over_baseline\\tgc_percent\\tmito\\n' > ${sample}.mitochondrion.candidates.tsv
    printf 'contig\\tdecision\\n' > ${sample}.mitochondrion.validation.tsv
    : > ${sample}.mitochondrion.depth.txt
    : > ${sample}.mitochondrion.warnings.txt
    """
}

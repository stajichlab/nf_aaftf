process ORGANELLE {
    tag   { sample }
    label 'aaftf_lite'
    publishDir "${params.outdir}/mitochondrion", mode: 'copy',
        pattern: '*.mitochondrion.{fasta.gz,candidates.tsv,validation.tsv,depth.txt,warnings.txt}'
    publishDir "${params.outdir}/plastid", mode: 'copy',
        pattern: '*.plastid.{fasta.gz,validation.tsv}'

    input:
    tuple val(sample), path(asm), path(filtered_1), path(filtered_2)

    output:
    tuple val(sample), path("${sample}.nuclear.fasta"), path("${sample}.mitochondrion.depth.txt"), emit: nuclear
    tuple val(sample), path("${sample}.mitochondrion.warnings.txt"), emit: warnings
    path "${sample}.mitochondrion.fasta.gz"
    path "${sample}.mitochondrion.candidates.tsv"
    path "${sample}.mitochondrion.validation.tsv"
    path "${sample}.plastid.fasta.gz"
    path "${sample}.plastid.validation.tsv"

    script:
    def require_hit = (params.getOrDefault('mitochondrion_require_hit', true) as String).toBoolean()
    def do_mito     = !(params.getOrDefault('skip_mitochondrion', false) as String).toBoolean()
    def do_plast    = !(params.getOrDefault('skip_plastid', true) as String).toBoolean()
    def veto_cov    = params.getOrDefault('rdna_veto_cov', 0.2)
    """
    set -euo pipefail
    trap 'echo "ORGANELLE: ERROR: command failed at line \$LINENO of .command.sh (exit \$?)" >&2' ERR
    # (the aaftf profile sets shell=/bin/bash without -ue, so fail fast explicitly)
    # Fail early with a clear message if this AAFTF image predates the depth MITO flag.
    depth_help=\$(AAFTF depth -h 2>&1 || true)
    if ! grep -q -- '--mito_fold' <<< "\$depth_help"; then
        echo "ORGANELLE: ERROR: this AAFTF image has no 'depth --mito_fold' option; use an image built from AAFTF main >= 0904010 (set --aaftf_sif)" >&2
        exit 1
    fi
    DB=${params.organelle_db}
    VALIDATE=${projectDir}/bin/organelle_validate.py
    # fasta_of IDS FASTA > subset
    fasta_of() { awk -v ids="\$1" 'BEGIN{while((getline l < ids)>0) k[l]=1} /^>/{id=substr(\$1,2); p=(id in k)} p' "\$2"; }
    # install_hint KIND
    install_hint() { echo "module load ncbi-blast/2.16.0+ && ${projectDir}/bin/download_organelle_db.sh -o \$DB \$1"; }

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
        echo "ORGANELLE: ERROR: AAFTF depth wrote no \$TSV; the container needs the MITO flag (AAFTF main >= 0904010)" >&2
        exit 1
    fi
    mv "\$TSV" ${sample}.mitochondrion.candidates.tsv

    # 2. Candidate pool. Mitochondrial candidates are the contigs the depth/GC (MITO) rule flagged. Plastid
    #    candidates (only when the plastid step is on) are contigs at >= plastid_min_fold x the nuclear baseline
    #    and >= mito_min_len bp. Only the pool is searched, never the whole assembly.
    : > pool.ids
    if [ "${do_mito}" = "true" ]; then
        awk -F'\\t' 'NR>1 && \$6=="yes" {print \$1}' ${sample}.mitochondrion.candidates.tsv >> pool.ids
    fi
    if [ "${do_plast}" = "true" ]; then
        awk -F'\\t' -v f=${params.plastid_min_fold} -v l=${params.mito_min_len} 'NR>1 && \$2>=l && \$4>=f {print \$1}' ${sample}.mitochondrion.candidates.tsv >> pool.ids
    fi
    sort -u pool.ids > pool.sorted && mv pool.sorted pool.ids
    fasta_of pool.ids ${asm} > pool.fasta

    # 3. Score the pool by whole-contig coverage (union of alignments / contig length, so a NUMT or NUPT
    #    insertion in a nuclear contig does not count).
    MITO_ARGS=""; PLAST_ARGS=""; RDNA_ARGS=""
    if [ "${do_mito}" = "true" ] && [ "${require_hit}" = "true" ]; then
        if [ ! -s "\$DB/mitochondrion/mitochondrion.nhr" ]; then
            echo "ORGANELLE: ERROR: no mitochondrion BLAST DB at \$DB/mitochondrion; run: \$(install_hint mitochondrion)  (or set --mitochondrion_require_hit false)" >&2
            exit 1
        fi
        \$VALIDATE --candidates pool.fasta --db \$DB/mitochondrion/mitochondrion --min-cov 2 \\
            --min-pident ${params.organelle_min_pident} --cpus ${task.cpus} --out mito.tsv --validated /dev/null
        MITO_ARGS="--mito mito.tsv"
    elif [ "${do_mito}" = "true" ]; then
        MITO_ARGS="--assume-mito"
    fi
    # The plastid database is used when it is installed: for plastid calls, and as a check that stops a
    # chloroplast contig from being removed as mitochondrial. Without it, mitochondrial calls work as before.
    if [ -s "\$DB/plastid/plastid.nhr" ]; then
        \$VALIDATE --candidates pool.fasta --db \$DB/plastid/plastid --min-cov 2 \\
            --min-pident ${params.organelle_min_pident} --cpus ${task.cpus} --out plast.tsv --validated /dev/null
        PLAST_ARGS="--plast plast.tsv"
    elif [ "${do_plast}" = "true" ]; then
        echo "ORGANELLE: ERROR: no plastid BLAST DB at \$DB/plastid; run: \$(install_hint plastid)  (or set --skip_plastid true)" >&2
        exit 1
    else
        echo "ORGANELLE: note: no plastid database at \$DB/plastid, so mitochondrial calls are not checked against plastid genomes" >&2
    fi
    # rDNA veto: nuclear rDNA is high copy and can look like an organelle. It also labels the depth
    # OUTLIER contigs in the warnings (rDNA-like or unexplained).
    awk '/\\*\\* OUTLIER/ && !/MITO/ {print \$1}' ${sample}.mitochondrion.depth.txt > outlier.ids
    cat pool.ids outlier.ids | sort -u > rdna_query.ids
    : > rdna.tsv
    if [ -s rdna_query.ids ] && [ "${veto_cov}" != "0" ]; then
        if [ ! -s "\$DB/rdna/rdna.nal" ]; then
            echo "ORGANELLE: ERROR: no rDNA BLAST DB at \$DB/rdna; run: \$(install_hint rdna)  (or set --rdna_veto_cov 0)" >&2
            exit 1
        fi
        fasta_of rdna_query.ids ${asm} > rdna_query.fasta
        \$VALIDATE --candidates rdna_query.fasta --db \$DB/rdna/rdna --min-cov 2 --min-pident 80 \\
            --cpus ${task.cpus} --out rdna.tsv --validated /dev/null
    fi
    if [ -s rdna.tsv ]; then RDNA_ARGS="--rdna rdna.tsv"; fi

    # 4. Decide (bin/organelle_classify.py, unit-tested): rDNA veto, the better database wins, a tie goes to
    #    mitochondrion only if the MITO rule fired, each class needs its own minimum coverage.
    REMOVE_ARGS=""
    if [ "${do_mito}" = "true" ]; then REMOVE_ARGS="\$REMOVE_ARGS --remove-mito"; fi
    if [ "${do_plast}" = "true" ]; then REMOVE_ARGS="\$REMOVE_ARGS --remove-plastid"; fi
    ${projectDir}/bin/organelle_classify.py --candidates ${sample}.mitochondrion.candidates.tsv --pool pool.ids \\
        \$MITO_ARGS \$PLAST_ARGS \$RDNA_ARGS \$REMOVE_ARGS \\
        --mito-min-cov ${params.mito_min_cov} --plastid-min-cov ${params.plastid_min_cov} --rdna-veto-cov ${veto_cov} \\
        --out ${sample}.mitochondrion.validation.tsv --remove-mito-ids mito.rm --remove-plastid-ids plast.rm \\
        --warnings ${sample}.mitochondrion.warnings.txt
    cat mito.rm plast.rm > removed.ids
    # plastid table: every candidate the classifier assigned to the plastid class (removed or only reported)
    awk -F'\\t' 'NR==1 || \$13=="plastid"' ${sample}.mitochondrion.validation.tsv > ${sample}.plastid.validation.tsv

    # 5. Split three ways: mitochondrion, plastid, and everything else continues.
    awk -v ids=mito.rm 'BEGIN{while((getline l < ids)>0) k[l]=1} /^>/{id=substr(\$1,2); p=(id in k)} p' ${asm} > ${sample}.mitochondrion.fasta
    awk -v ids=plast.rm 'BEGIN{while((getline l < ids)>0) k[l]=1} /^>/{id=substr(\$1,2); p=(id in k)} p' ${asm} > ${sample}.plastid.fasta
    awk -v ids=removed.ids 'BEGIN{while((getline l < ids)>0) k[l]=1} /^>/{id=substr(\$1,2); p=!(id in k)} p' ${asm} > ${sample}.nuclear.fasta
    n_in=\$(grep -c '^>' ${asm} || true); n_nuc=\$(grep -c '^>' ${sample}.nuclear.fasta || true)
    n_mt=\$(grep -c '^>' ${sample}.mitochondrion.fasta || true); n_pl=\$(grep -c '^>' ${sample}.plastid.fasta || true)
    if [ "\$((n_nuc + n_mt + n_pl))" -ne "\$n_in" ] || [ "\$n_nuc" -eq 0 ]; then
        echo "ORGANELLE: ERROR: split mismatch for ${sample}: in=\$n_in nuclear=\$n_nuc mitochondrion=\$n_mt plastid=\$n_pl" >&2
        exit 1
    fi
    echo "ORGANELLE: ${sample}: \$n_in contigs -> \$n_nuc kept, \$n_mt mitochondrial and \$n_pl plastid removed (\$(wc -l < pool.ids) candidates)"
    gzip -c ${sample}.mitochondrion.fasta > ${sample}.mitochondrion.fasta.gz
    gzip -c ${sample}.plastid.fasta > ${sample}.plastid.fasta.gz

    # 6. Warnings (nothing removed for these): candidates that were kept (written above), and depth OUTLIER
    #    contigs that were not removed, labeled by rDNA coverage.
    awk -v OFS='\\t' -v veto=${veto_cov} '
        FILENAME=="removed.ids" { rm[\$1]=1; next }
        FILENAME=="rdna.tsv" { if (FNR>1) rd[\$1]=\$4; next }
        /\\*\\* OUTLIER/ && !(\$1 in rm) { r = (\$1 in rd) ? rd[\$1] : 0
            print "OUTLIER", \$1, "len=" \$2, "depth=" \$3, "GC=" \$4, "rdna_cov=" r, (r>=veto ? "rDNA-like" : "unexplained") }' removed.ids rdna.tsv ${sample}.mitochondrion.depth.txt >> ${sample}.mitochondrion.warnings.txt
    if [ -s ${sample}.mitochondrion.warnings.txt ]; then
        echo "ORGANELLE: WARNING: ${sample}: contigs kept that need review:" >&2
        cat ${sample}.mitochondrion.warnings.txt >&2
    fi
    ${projectDir}/bin/record_versions.sh "${params.outdir}/sample_versions/${sample}/versions" "${task.process}" "${sample}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'minimap2=minimap2 --version' \\
        'samtools=samtools --version | head -1' \\
        'mosdepth=mosdepth --version' \\
        'blastn=blastn -version | head -1' \\
        "mitochondrion_db=cat ${params.organelle_db}/mitochondrion/RELEASE_NUMBER 2>/dev/null || echo unknown" \\
        "plastid_db=cat ${params.organelle_db}/plastid/RELEASE_NUMBER 2>/dev/null || echo not_installed"
    """

    stub:
    """
    cp ${asm} ${sample}.nuclear.fasta
    : > ${sample}.mitochondrion.fasta
    gzip -c ${sample}.mitochondrion.fasta > ${sample}.mitochondrion.fasta.gz
    : > ${sample}.plastid.fasta
    gzip -c ${sample}.plastid.fasta > ${sample}.plastid.fasta.gz
    printf 'contig\\tlength\\tmean_depth\\tfold_over_baseline\\tgc_percent\\tmito\\n' > ${sample}.mitochondrion.candidates.tsv
    printf 'contig\\tdecision\\n' > ${sample}.mitochondrion.validation.tsv
    printf 'contig\\tdecision\\n' > ${sample}.plastid.validation.tsv
    : > ${sample}.mitochondrion.depth.txt
    : > ${sample}.mitochondrion.warnings.txt
    """
}

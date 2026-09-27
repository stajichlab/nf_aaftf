#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * nf_aaftf — AAFTF genome assembly + cleanup for short-read (Illumina) genomes.
 *
 * Flow:  TRIM -> FILTER -> ASSEMBLE -> [vector_screen] -> [contamination_screen]
 *        -> RMDUP -> POLISH -> SORT -> {COMPRESS, ASSESS, [optional DEPTH]}
 *
 * The final assembly fasta is published bgzip-compressed (results/sort/*.sorted.fasta.gz)
 * so it stays indexable via samtools/bcftools while saving disk space.
 *
 * Screening stages (both optional, both complementary):
 *   vector_screen:       vecscreen (BLASTN, default) OR fcs_screen (NCBI FCS adaptor)
 *   contamination_screen: fcs_gx (NCBI FCS-GX purge) AND/OR sourpurge (sourmash purge)
 *
 * All processes run inside the AAFTF container (singularity SIF or Docker image,
 * selected via params.container_engine). Focus: short-read assembly with SPAdes
 * + POLCA polishing. See conf/profile_aaftf.config for tunables and README.md for
 * design notes (JGI BBMap best-practices, alternative assemblers).
 */

params.samples = "${launchDir}/samples.csv"
params.indir   = "${launchDir}/input"
params.outdir  = "${launchDir}/results"
params.n_test  = 0

// ── Include process modules ─────────────────────────────────────────────────
include { AAFTF_TRIM } from './modules/aaftf/TRIM'
include { FILTER     } from './modules/aaftf/FILTER'
include { ASSEMBLE   } from './modules/aaftf/ASSEMBLE'
include { VECSCREEN  } from './modules/aaftf/VECSCREEN'
include { FCS_SCREEN } from './modules/aaftf/FCS_SCREEN'
include { CONTAM_CLEAN } from './modules/aaftf/CONTAM_CLEAN'
include { CONTAM_CLEAN_BATCH } from './modules/aaftf/CONTAM_CLEAN_BATCH'
include { SOURPURGE  } from './modules/aaftf/SOURPURGE'
include { RMDUP      } from './modules/aaftf/RMDUP'
include { POLISH     } from './modules/aaftf/POLISH'
include { SORT       } from './modules/aaftf/SORT'
include { COMPRESS   } from './modules/aaftf/COMPRESS'
include { ASSESS     } from './modules/aaftf/ASSESS'
include { DEPTH      } from './modules/aaftf/DEPTH'
// Published intermediate FASTA (asm, vecscreen/fcs_screen, contam_clean) are
// gzip-compressed to save disk space, but the durable-resume paths below
// read those published files back in as real pipeline data and AAFTF's own
// FASTA I/O can't read gzip -- so each reuse path decompresses via its own
// GUNZIP invocation (one process can't be included/called under 3 names
// otherwise, since Nextflow requires a distinct alias per call site).
include { GUNZIP as GUNZIP_ASM_REUSE } from './modules/aaftf/GUNZIP'
include { GUNZIP as GUNZIP_VEC_REUSE } from './modules/aaftf/GUNZIP'
include { GUNZIP as GUNZIP_CONTAM    } from './modules/aaftf/GUNZIP'

// ── Durable, file-existence-based stage-reuse helpers ───────────────────────
// Nextflow's own -resume cache is fragile: a corrupted/reset cache database
// (confirmed to happen in practice, 2026-09-26 -- an unrelated `rm -rf
// .nextflow*` in the same launchDir silently reset this pipeline's session
// cache) forces every stage to recompute even though the real outputs are
// still safely on disk. These checks are independent of -resume/cache
// health: a stage is skipped whenever ITS OWN final published output already
// exists in outdir. Declared as top-level functions (not closures assigned
// to `def` locals) because Nextflow's DSL2 script compiler does not resolve
// one local closure calling another from inside an operator closure
// (e.g. `.branch{}`) -- confirmed by compile failure when tried that way.
def allExist(List fs) {
    fs.every { it.exists() && it.size() > 0 }
}

def vecOutDir() {
    (params.getOrDefault('vector_screen_method', 'vecscreen') == 'fcs_screen') ? 'fcs_screen' : 'vecscreen'
}

def vecExt() {
    (params.getOrDefault('vector_screen_method', 'vecscreen') == 'fcs_screen') ? 'fcs_screen.fasta' : 'vecscreen.fasta'
}

def filterDone(String s) {
    allExist([
        file("${params.outdir}/filter/${s}_filtered_1.fastq.gz"),
        file("${params.outdir}/filter/${s}_filtered_2.fastq.gz"),
        file("${params.outdir}/filter/${s}_filtered_U.fastq.gz")
    ])
}

def vecDone(String s, boolean skipVecscreen) {
    if (skipVecscreen) { return false }
    // Published copy is gzip-compressed (see VECSCREEN/FCS_SCREEN).
    def f = file("${params.outdir}/${vecOutDir()}/${s}.${vecExt()}.gz")
    f.exists() && f.size() > 0
}

// ── Workflow ────────────────────────────────────────────────────────────────
workflow {
    main:

    // Boolean CLI params arrive as strings ("true"/"false"), so `--skip_fcsgx
    // false` would otherwise be truthy in Groovy. Coerce explicitly here.
    def skip_vecscreen = (params.getOrDefault('skip_vecscreen', false) as String).toBoolean()
    def skip_fcsgx     = (params.getOrDefault('skip_fcsgx', true)     as String).toBoolean()
    def skip_sourpurge = (params.getOrDefault('skip_sourpurge', true) as String).toBoolean()
    def run_depth      = (params.getOrDefault('run_depth', true)      as String).toBoolean()
    def vec_method     = params.getOrDefault('vector_screen_method', 'vecscreen')

    // --reuse_asm V11,V20: take these samples' assemblies from
    // ${outdir}/asm/<sample>.<assembler>.fasta instead of running ASSEMBLE.
    // Use it when a finished assembly was marked failed (e.g. a job moved
    // between partitions) and -resume would otherwise re-assemble it. Only
    // the listed samples switch source, so other samples keep their
    // downstream cache (the cache hash includes the input file path).
    def reuse_param = params.getOrDefault('reuse_asm', '')
    if (reuse_param instanceof Boolean || (reuse_param as String) == 'true') {
        error("--reuse_asm needs a comma-separated sample list, e.g. --reuse_asm V11,V20")
    }
    def reuse_set = ((reuse_param ?: '') as String).tokenize(',')*.trim().findAll { it } as Set
    reuse_set.each { s ->
        // Published copy is gzip-compressed (see ASSEMBLE).
        def f = file("${params.outdir}/asm/${s}.${params.assembler}.fasta.gz")
        if (!f.exists()) {
            error("--reuse_asm: no assembly for ${s} at ${f}")
        }
    }
    if (reuse_set) {
        log.info "Reusing existing assemblies (ASSEMBLE skipped): ${reuse_set.sort().join(', ')}"
    }

    // ── Durable, file-existence-based stage reuse ────────────────────────
    // See allExist/vecOutDir/vecExt/filterDone/vecDone above for the
    // mechanism. Two independent axes, because filtered reads are needed
    // downstream (POLISH/DEPTH) regardless of whether the assembly/vecscreen
    // side is reused:
    //   Axis 1 (reads):    FILTER done?  -> skip TRIM+FILTER, reuse filtered reads.
    //   Axis 2 (assembly): VECSCREEN/FCS_SCREEN done? -> skip ASSEMBLE+screen entirely.

    // Synchronous pre-scan (samplesheet parsed directly, not via the
    // reactive Channel below) purely so the two reuse axes get the same
    // up-front log visibility as --reuse_asm above.
    def all_sample_ids = file(params.samples).splitCsv(header: true, sep: ',').collect { it.sample.trim() }
    def filter_reuse_ids = all_sample_ids.findAll { filterDone(it) }
    def vec_reuse_ids    = all_sample_ids.findAll { vecDone(it, skip_vecscreen) }
    if (filter_reuse_ids) {
        log.info "Reusing existing filtered reads (TRIM+FILTER skipped): ${filter_reuse_ids.sort().join(', ')}"
    }
    if (vec_reuse_ids) {
        log.info "Reusing existing screened assemblies (ASSEMBLE+${vecOutDir().toUpperCase()} skipped): ${vec_reuse_ids.sort().join(', ')}"
    }

    // Sample sheet: sample,read_1,read_2 [,taxid]. The optional taxid column
    // feeds the optional FCS-GX / sourpurge steps (NCBI taxonomy id, e.g.
    // 4751 Fungi / 4890 Ascomycota); falls back to params.fcs_taxid when absent.
    Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def tax = row.taxid ? row.taxid.trim() : params.fcs_taxid
            tuple(row.sample.trim(),
                  file("${params.indir}/${row.read_1}", checkIfExists: true),
                  file("${params.indir}/${row.read_2}", checkIfExists: true),
                  tax)
        }
        .ifEmpty { error("No samples found in ${params.samples}") }
        .take( (params.n_test as int) > 0 ? (params.n_test as int) : Integer.MAX_VALUE )
        .set { samples_ch }

    // ── Stage 1+2: QC trim + contaminant read filtering ──────────────
    // Axis 1 (reads): skip TRIM+FILTER entirely for samples whose filtered
    // reads are already published (durable, file-existence-based -- see
    // filterDone() above). Filtered reads are needed downstream (POLISH,
    // DEPTH) regardless of what happens on the assembly/vecscreen axis, so
    // this check is independent of reuse_set / vecDone.
    def ch_samples_split = samples_ch.branch { s, r1, r2, t ->
        filter_reuse: filterDone(s)
        fresh:        true
    }
    AAFTF_TRIM(ch_samples_split.fresh.map { s, r1, r2, t -> tuple(s, r1, r2) })
    FILTER(AAFTF_TRIM.out.trimmed)
    def ch_filtered = FILTER.out.filtered.mix(
        ch_samples_split.filter_reuse.map { s, r1, r2, t ->
            tuple(s,
                  file("${params.outdir}/filter/${s}_filtered_1.fastq.gz"),
                  file("${params.outdir}/filter/${s}_filtered_2.fastq.gz"),
                  file("${params.outdir}/filter/${s}_filtered_U.fastq.gz"))
        }
    )

    // ── Stage 3+4: assembly (SPAdes) + vector/primer screening ───────
    // Axis 2 (assembly): skip ASSEMBLE *and* VECSCREEN/FCS_SCREEN entirely
    // for samples whose final screened assembly is already published
    // (vecDone() above) -- this is the common case once a sample has fully
    // cleared this stage in a prior run. --reuse_asm (explicit sample list)
    // still applies for samples that need a fresh screen off an existing,
    // not-yet-screened assembly.
    def ch_vec_split = ch_filtered.branch { s, f1, f2, fu ->
        vec_reuse: vecDone(s, skip_vecscreen)
        process:   true
    }
    def ch_filt = ch_vec_split.process.branch { s, f1, f2, fu ->
        asm_reuse: s in reuse_set
        assemble:  true
    }
    ASSEMBLE(ch_filt.assemble)
    GUNZIP_ASM_REUSE(
        ch_filt.asm_reuse.map { s, f1, f2, fu ->
            tuple(s, file("${params.outdir}/asm/${s}.${params.assembler}.fasta.gz"))
        }
    )
    def ch_asm = ASSEMBLE.out.assembly.mix(GUNZIP_ASM_REUSE.out.plain)

    //   vecscreen (BLASTN, default) OR fcs_screen (NCBI FCS adaptor)
    def ch_vec_fresh
    if (skip_vecscreen) {
        ch_vec_fresh = ch_asm
    } else if (vec_method == 'fcs_screen') {
        FCS_SCREEN(ch_asm)
        ch_vec_fresh = FCS_SCREEN.out.screened
    } else {
        VECSCREEN(ch_asm)
        ch_vec_fresh = VECSCREEN.out.vecscreen
    }
    GUNZIP_VEC_REUSE(
        ch_vec_split.vec_reuse.map { s, f1, f2, fu ->
            tuple(s, file("${params.outdir}/${vecOutDir()}/${s}.${vecExt()}.gz"))
        }
    )
    def ch_vec = ch_vec_fresh.mix(GUNZIP_VEC_REUSE.out.plain)

    // ── Stage 4b: contamination screening ───────────────────────────
    //   fcs_gx (NCBI FCS-GX purge) AND/OR sourpurge (sourmash purge)
    // Both can run; fcs_gx first then sourpurge, or either alone.
    def ch_purged = ch_vec
    if (!skip_fcsgx) {
        // The FCS-GX database (~465 GB) is far too expensive to rsync-stage
        // once per genome, so above contam_clean_batch_size=0 we batch
        // genomes into single SLURM jobs that stage the DB once and clean
        // every genome in the batch (~15-30 min staging + ~1-2 min/genome).
        // Set contam_clean_batch_size = 0 to fall back to the original
        // one-job-per-genome CONTAM_CLEAN process.
        int contam_clean_batch_size = params.getOrDefault('contam_clean_batch_size', 100) as int
        def items_ch = ch_vec.join(samples_ch.map { s, r1, r2, t -> tuple(s, t) })
            .map { s, f, t -> tuple(s, f.name, t, f) }

        if (contam_clean_batch_size > 0) {
            // Skip samples a prior attempt already cleaned, so re-launching
            // the pipeline never repays the DB-staging cost for them.
            // Published copy is gzip-compressed (see CONTAM_CLEAN_BATCH).
            def items_to_clean = items_ch.filter { s, fname, t, f ->
                !file("${params.outdir}/contam_clean/${s}.contam_clean.fasta.gz").exists()
            }
            def batches = items_to_clean.collate(contam_clean_batch_size)
                .map { batch ->
                    tuple(batch.collect { s, fname, t, f -> [s, fname, t] },
                          batch.collect { s, fname, t, f -> f })
                }
            CONTAM_CLEAN_BATCH(batches)
            def clean_done_ch = CONTAM_CLEAN_BATCH.out.manifest.collect().ifEmpty([])

            // The cleaned assembly always lands at the fixed path below
            // (whether cleaned just now or in a prior attempt), so rebuild
            // the per-sample channel from that convention rather than from
            // CONTAM_CLEAN_BATCH's own (batch-shaped) output. It's gzipped
            // on disk (space), so decompress it back to plain FASTA before
            // handing it to RMDUP -- every sample goes through this, not
            // just resumed ones, since CONTAM_CLEAN_BATCH always writes .gz.
            def ch_contam_gz = ch_vec
                .map { s, f -> tuple(s, file("${params.outdir}/contam_clean/${s}.contam_clean.fasta.gz")) }
                .combine(clean_done_ch)
                .map { it[0..1] }
                .filter { s, f ->
                    if (!f.exists()) {
                        log.warn "CONTAM_CLEAN_BATCH: no cleaned assembly for ${s}; skipping downstream"
                        return false
                    }
                    return true
                }
            GUNZIP_CONTAM(ch_contam_gz)
            ch_purged = GUNZIP_CONTAM.out.plain
        } else {
            CONTAM_CLEAN(items_ch.map { s, fname, t, f -> tuple(s, f, t) })
            ch_purged = CONTAM_CLEAN.out.clean
        }
    }
    if (!skip_sourpurge) {
        SOURPURGE(ch_purged.join(samples_ch.map { s, r1, r2, t -> tuple(s, t) }))
        ch_purged = SOURPURGE.out.clean
    }

    // ── Stage 5: finishing ──────────────────────────────────────────
    // NOTE: must use ch_filtered (fresh + filter_reuse merged), not
    // FILTER.out.filtered directly -- that channel now only carries the
    // freshly-filtered samples, since filter_reuse samples bypass FILTER
    // entirely above.
    RMDUP(ch_purged)
    POLISH(RMDUP.out.rmdup.join(ch_filtered.map { s, f1, f2, fu -> tuple(s, f1, f2) }))
    SORT(POLISH.out.polished)
    COMPRESS(SORT.out.sorted)
    ASSESS(SORT.out.sorted)

    if (run_depth) {
        DEPTH(SORT.out.sorted.join(ch_filtered.map { s, f1, f2, fu -> tuple(s, f1, f2) }))
    }
}

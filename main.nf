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

    // ── Stage 1: QC trim (fastp dedup+merge+cuttail, then bbduk) ───
    AAFTF_TRIM(samples_ch.map { s, r1, r2, t -> tuple(s, r1, r2) })

    // ── Stage 2: contaminant read filtering (phiX + vector DBs) ─────
    FILTER(AAFTF_TRIM.out.trimmed)

    // ── Stage 3: assembly (SPAdes) ──────────────────────────────────
    ASSEMBLE(FILTER.out.filtered)

    // ── Stage 4: vector / primer screening ──────────────────────────
    //   vecscreen (BLASTN, default) OR fcs_screen (NCBI FCS adaptor)
    def ch_vec
    if (skip_vecscreen) {
        ch_vec = ASSEMBLE.out.assembly
    } else if (vec_method == 'fcs_screen') {
        FCS_SCREEN(ASSEMBLE.out.assembly)
        ch_vec = FCS_SCREEN.out.screened
    } else {
        VECSCREEN(ASSEMBLE.out.assembly)
        ch_vec = VECSCREEN.out.vecscreen
    }

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
            def items_to_clean = items_ch.filter { s, fname, t, f ->
                !file("${params.outdir}/contam_clean/${s}.contam_clean.fasta").exists()
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
            // CONTAM_CLEAN_BATCH's own (batch-shaped) output.
            ch_purged = ch_vec
                .map { s, f -> tuple(s, file("${params.outdir}/contam_clean/${s}.contam_clean.fasta")) }
                .combine(clean_done_ch)
                .map { it[0..1] }
                .filter { s, f ->
                    if (!f.exists()) {
                        log.warn "CONTAM_CLEAN_BATCH: no cleaned assembly for ${s}; skipping downstream"
                        return false
                    }
                    return true
                }
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
    RMDUP(ch_purged)
    POLISH(RMDUP.out.rmdup.join(FILTER.out.filtered.map { s, f1, f2, fu -> tuple(s, f1, f2) }))
    SORT(POLISH.out.polished)
    COMPRESS(SORT.out.sorted)
    ASSESS(SORT.out.sorted)

    if (run_depth) {
        DEPTH(SORT.out.sorted.join(FILTER.out.filtered.map { s, f1, f2, fu -> tuple(s, f1, f2) }))
    }
}

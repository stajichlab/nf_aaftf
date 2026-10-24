# nf_aaftf

AAFTF genome **assembly and cleanup** for short-read (Illumina) genomes, as a
Nextflow DSL2 pipeline that runs **exclusively inside the AAFTF container**
(Singularity SIF or the same image via Docker). SPAdes assembly + POLCA
polishing, with `fcs-gx`, vector screening, and depth estimation as optional
steps.

This is the Nextflow port of the shell pipeline in
`Desert_microbes/Culture_genomes/pipeline/AAFTF/` (`01_AAFTF.sh`,
`02_clean_fcsgx.sh`, `03_AAFTF_finish.sh`).

## Flow

```
samples.csv
 │
 ├─ TRIM ───────────  pass 1: fastp (default: dedup+merge+cuttail) │ trimmomatic
 │    └─ pass 2: bbduk (fixed second adaptor/quality pass)
 │
 ├─ FILTER ──────────  bbduk (default) │ bowtie2 │ bwa │ minimap2
 │                     vs phiX + contam/vector DBs
 │
 ├─ ASSEMBLE ────────  spades (default) │ megahit │ unicycler
 │
 ├─ VECSCREEN ───────  vector/primer screening (optional, ON by default)
 │    ├─ fcs_screen (default) — NCBI FCS adaptor screening
 │    └─ vecscreen            — BLASTN vs UniVec + contaminant DBs
 │
 ├─ CONTAM_CLEAN ────  contamination screening (optional, OFF by default)
 │    ├─ fcs_gx (default)     — NCBI FCS-GX taxonomy-based purge
 │    └─ sourpurge (optional) — AAFTF taxpurge: sourmash taxonomy screen
 │
 ├─ RMDUP ───────────  remove duplicate contigs
 │
 ├─ POLISH ──────────  polca (default) │ pypolca │ polypolish │ nextpolish2 │ racon
 │
 ├─ SORT ────────────  by length, drop < min_contig_len, rename headers
 │
 ├─ ASSESS ──────────  assembly stats  (always)
 │
 └─ DEPTH ───────────  read depth back-mapping (optional, ON by default)
```

**Stage tool choices** — each maps to an AAFTF `--method` flag:

| Stage | Param | Options (bold = default) | Notes |
|-------|-------|--------------------------|-------|
| TRIM | `trim_method1` | **fastp** \| trimmomatic | Primary QC + adaptor trim; fastp adds dedup+merge |
| TRIM | `trim_method2` | **bbduk** | Second adaptor pass (fixed; uses interleave workaround) |
| FILTER | `filter_aligner` | **bbduk** \| bowtie2 \| bwa \| minimap2 | Read filter aligner vs phiX + contam/vector DBs |
| ASSEMBLE | `assembler` | **spades** \| megahit \| unicycler | dipspades was removed in AAFTF 0.7 |
| ASSEMBLE | `assembler_args` | (string, default null) | Extra CLI args passed to the assembler |
| VECSCREEN | `vector_screen_method` | **fcs_screen** \| vecscreen | Vector/primer screening method |
| CONTAM_CLEAN | `skip_fcsgx` | **true** \| false | NCBI FCS-GX contamination purge |
| CONTAM_CLEAN | `skip_sourpurge` | **true** \| false | Sourmash-based contamination purge |
| POLISH | `polisher` | **polca** \| pypolca \| polypolish \| nextpolish2 \| racon | POLCA uses original filtered reads |

Outputs land under `results/<step>/` (publishDir), one file set per sample.

## Container

The pipeline is container-only. Set a single engine flag:

| `params.container_engine` | Container used |
|---|---|
| `singularity` (default) | cached SIF `AAFTF-main-ccd6c26.sif` (pulled from `ghcr.io/stajichlab/aaftf:main`, commit ccd6c26; pinned until a stable 0.7.0 release) |
| `docker`                 | `docker://ghcr.io/stajichlab/aaftf:latest` |

The AAFTF image carries its toolchain (`spades.py`, `megahit`, `bbduk.sh`, `fastp`, `minimap2`,
`bwa`, `bowtie2`, `samtools`, `mosdepth`, `polypolish`, `pypolca`, `aaftf_polca.sh`, `sourmash` ...;
run `AAFTF dependency` in the image for the full list). The pipeline binds the host reference
DB (at `/opt/aaftf_db`) and the NCBI taxdump into the container via `singularity.runOptions`.
**`taxonkit` is not in the AAFTF 0.7 image** (removed in commit `1f51deb`). The `TAXON_PHYLUM`
process runs it from its own container, `params.taxonkit_container` (default
`docker://quay.io/biocontainers/taxonkit:0.20.0--h9ee0642_1`). It gives the phylum taxid to
CONTAM_CLEAN and CONTAM_CLEAN_BATCH and the phylum name to SOURPURGE. It runs only when FCS-GX or
sourpurge is enabled.

> Note: `AAFTF --version` prints `0.0.0+unknown` for this dev (`latest`) build —
> cosmetic only; all subcommands and tools work (verified).

## Verify the container works (what I smoke-tested)

Checked inside `aaftf_latest.sif`:
- `AAFTF trim --method fastp --dedup --merge --cutright` → `_1P/_2P/_MG.fastq.gz` ✅
- `AAFTF filter --aligner bbduk` with `AAFTF_DB` bound to
  `/bigdata/stajichlab/shared/lib/AAFTF_DB` (the module default
  `/srv/projects/db/AAFTF_DB` is **missing phiX** — that's the earlier failure) ✅
- `AAFTF assemble` (spades) — failed only on a 10k-read subsample where all
  paired reads were filtered (data artifact). With 1.2M-read test data the real
  run assembled 5 contigs end-to-end through depth (see below).

## Quick start

```bash
# 1. (optional) validate the graph fast, no real compute:
nextflow run . -profile test -stub-run --n_test 2

# 2. real run (reads must be in input/, sample sheet at samples.csv)
sbatch run_aaftf.sh
#    -- OR run directly on a node:
nextflow run . -profile aaftf -resume
```

`samples.csv` columns: `sample,read_1,read_2[,taxid]`. The optional `taxid`
column is the NCBI taxonomy id of the organism (e.g. `4890` Ascomycota / `4751`
Fungi) and feeds the optional FCS-GX step; falls back to `params.fcs_taxid`.

## Configuration

Everything tunable lives in `conf/profile_aaftf.config`. The key knobs are the
**per-stage tool choices**, each of which maps to an AAFTF `--method` flag:

### Stage tool choices

| Param | Options (bold = default) | What it controls |
|-------|--------------------------|------------------|
| `trim_method1` | **fastp** \| trimmomatic | Primary QC/adaptor trim pass |
| `trim_method2` | **bbduk** | Second adaptor/quality pass (fixed) |
| `filter_aligner` | **bbduk** \| bowtie2 \| bwa \| minimap2 | Read filter aligner |
| `assembler` | **spades** \| megahit \| unicycler | Assembly method |
| `assembler_args` | string (default null) | Extra CLI args for the assembler |
| `polisher` | **polca** \| pypolca \| polypolish \| nextpolish2 \| racon | Polishing tool |

### Trimming parameters

When `trim_method1 = 'fastp'` (default): dedup + merge + 3' quality cut
(`--dedup --merge --cuttail`). When `trim_method1 = 'trimmomatic'`: uses the
`trimmomatic_*` knobs (adaptors, clip, leading/trailing/sliding window, quality).

| Param | Default | Description |
|-------|---------|-------------|
| `dedup` | true | fastp read deduplication |
| `minlen` | 75 | Minimum read length after trimming |
| `PHRED` | 20 | fastp 3' quality cutoff (cuttail) |
| `trim_memory` | 8 | GB for bbduk -Xmx (trim pass) |
| `trimmomatic_adaptors` | TruSeq3-PE.fa | Trimmomatic adaptor file |
| `trimmomatic_clip` | 2:30:10 | ILLUMINACLIP settings |
| `trimmomatic_leadingwindow` | 3 | LEADING quality threshold |
| `trimmomatic_trailingwindow` | 3 | TRAILING quality threshold |
| `trimmomatic_slidingwindow` | 4:15 | SLIDINGWINDOW:windowSize:requiredQuality |
| `trimmomatic_quality` | 15 | AVGQUAL threshold |

### Pipeline control flags

| Param | Default | Description |
|-------|---------|-------------|
| `skip_fcsgx` | true | FCS-GX purge (big DB download + highmem) |
| `skip_sourpurge` | true | Sourmash-based contamination purge |
| `skip_vecscreen` | false | Vector/primer screening on assembled contigs |
| `run_depth` | true | Read depth back-mapping |
| `reuse_asm` | '' | Comma-separated samples (e.g. `V11,V20`) whose assembly is taken from `<outdir>/asm/<sample>.<assembler>.fasta` instead of running ASSEMBLE. Use when a finished assembly was marked failed and `-resume` would re-assemble it. The run stops if a listed file is missing. |
| `fcs_taxid` | 4751 | Fallback NCBI taxonomy id for FCS-GX (4751 = Fungi) |

### Screening stages

Two complementary screening categories run after assembly:

**1. Vector / primer screening** (`skip_vecscreen` to disable):

| Param | Options (bold = default) | Description |
|-------|--------------------------|-------------|
| `vector_screen_method` | **fcs_screen** \| vecscreen | `fcs_screen` = NCBI FCS adaptor screening; `vecscreen` = BLASTN vs UniVec + contaminant DBs |
| `fcs_screen_prok` | **false** \| true | fcs_screen mode: `--euk` (default) or `--prok` |

**2. Contamination screening — CONTAM_CLEAN stage** (enable with `skip_fcsgx` / `skip_sourpurge`):

| Param | Default | Description |
|-------|---------|-------------|
| `skip_fcsgx` | **true** | NCBI FCS-GX taxonomy-based purge. Resolves phylum taxid from sample's organism taxid via `taxonkit` |
| `skip_sourpurge` | **true** | Sourmash k-mer LCA purge. Resolves phylum name from taxid and purges contigs not matching expected taxonomy |

Both contamination screens can run simultaneously (FCS-GX first, then sourpurge).
FCS-GX needs a current gxdb path (`params.fcsgx_db`) and high memory.

### Assembly / polishing / finish

| Param | Default | Description |
|-------|---------|-------------|
| `spades_memory` | 64 | GB for SPAdes |
| `polisher_iterations` | 1 | AAFTF polish --iterations (polca ignores) |
| `min_contig_len` | 2000 | AAFTF sort/rmdup min contig length |

Boolean CLI flags are coerced (`--skip_fcsgx false` works — a string `"false"`
would otherwise be truthy in Groovy).

## Design notes & best practices

### Read QC / filtering order (logic)
- **Trim first, filter second.** Trimming removes adapters + low-quality bases
  *before* screening reads against reference DBs, so spurious matches from
  adapter-contaminated reads don't cause false contamination calls.
- Two trimming passes mirror JGI-style BBDuk-first QC while staying inside
  AAFTF: `fastp` (dedup + merge + 3' quality trim) then `bbduk` for a second
  adapter/quality pass. Both use AAFTF's shipped default parameters.
  - **Confirmed BBMap 39.80 caveat** (in `AAFTF.sif`): `bbduk` crashes its
    `PairStreamer` with `List size mismatch` on variable-length **paired** reads
    (the fastp pass-1 output), silently dropping all reads if fed as two files.
    The TRIM module therefore feeds bbduk an **interleaved single file**
    (`shuffle.sh` → bbduk `interleaved=true` → `reformat.sh` back to pairs),
    which processes 100% of reads. Verified: 732,780 reads in → 727,924 kept.
- Reads are filtered against `phiX + contam/vector` DBs with the same aligner
  the user picks (`bbduk` default — the fast option; `bowtie2`/`bwa` are more
  sensitive but slower).
- Merged reads (from fastp `--merge`) are carried through the filter as the
  single-end `_U` set and given to SPAdes via `--merged`, so overlap-merged
  reads are not lost.

### Polishing / finishing
- POLCA is the default polisher (robust to indel error, uses original filtered
  reads), matching the source framework. Pilon/NextPolish/Racon are available.
- `sort` drops contigs `< min_contig_len` (default 2000) and renames headers to
  the sample name; `assess` reports N50 etc.; `depth` maps reads back with
  minimap2 and flags possible contaminant/organellar contigs by depth.

### Screening tools

AAFTF provides four screening subcommands in two complementary categories:

**Vector / primer screening** (run one):
- `fcs_screen` — NCBI FCS adaptor screening (default)
- `vecscreen` — BLASTN against UniVec + contaminant DBs (alternative)

Select with `vector_screen_method`. Disable entirely with `skip_vecscreen=true`.

**Contamination screening** (run either, both, or neither):
- `fcs_gx_purge` — NCBI FCS-GX taxonomy-based purge. Resolves the phylum taxid
  at runtime from the sample's organism taxid via `taxonkit`, then stages the
  FCS-GX database into node-local scratch before purging. Needs a current gxdb
  path (`params.fcsgx_db`) and high memory; off by default.
- `sourpurge` — runs `AAFTF taxpurge` (AAFTF >= 0.7 split `sourpurge` into `taxpurge` and
  `covpurge`; only the taxonomy screen is run here). Resolves the phylum name from the
  organism taxid and purges contigs whose taxonomy doesn't match the expected
  phylum. Needs the sourmash databases in `AAFTF_DB` (`AAFTF database taxpurge`, about
  1 GB for the default `-lite` databases); off by default.

Enable with `skip_fcsgx=false` and/or `skip_sourpurge=false`. Both can run in
sequence (FCS-GX → sourpurge).

### Alternative assemblers
- `assembler` is wired through to `AAFTF assemble --method`. Supported values:
  `spades` (default), `megahit`, `unicycler`. (`dipspades` was removed in AAFTF 0.7;
  the pipeline stops with a clear message if it is requested.)
- SPAdes, MEGAHIT and Unicycler ship in the AAFTF >= 0.7 image (`AAFTF dependency` lists them).
- Extra assembler-specific args can be passed via `assembler_args` (e.g.
  `--assembler_args '--kmer 21,33,55'` for spades/megahit).

## Layout

```
main.nf                  workflow (channel wiring + optional branches)
nextflow.config          manifest + SLURM executor + process defaults + profiles
conf/profile_aaftf.config  pipeline params + resources + workDir + trace/report
conf/test.config         -stub-run / tiny-profile validation
modules/aaftf/<STEP>/main.nf   one process per AAFTF step
run_aaftf.sh             sbatch launcher for the Nextflow head process
samples.csv              run sample sheet (reads in input/)
tests/data/              sample sheet; input/ FASTQs are too large for git (fetch separately)
```

## Validated end-to-end (2026-08-16)

Ran the full pipeline on the cluster with the real test data (`M40`,
`--n_test 1`, singularity, `aaftf_latest.sif`): 9/9 tasks succeeded
TRIM → FILTER → ASSEMBLE (SPAdes) → VECSCREEN → RMDUP → POLISH (POLCA) → SORT →
ASSESS → DEPTH. This surfaced and fixed three real-run issues:
partition time limits (2h cap on `short`), the BBMap `PairStreamer` bug
(interleave workaround), and a DEPTH publishDir glob. Outputs in
`tests/output_real/`.

## UCR HPCC notes

- Head process must run on a compute node → use `run_aaftf.sh` (never the login
  node for real runs).
- Partitions: `epyc` (default), `short` (trim/filter/lite), `highmem`
  (CONTAM_CLEAN / SPAdes escalation). Verify with `sinfo`.
- AAFTF image + reference DBs are lab-shared under
  `/bigdata/stajichlab/shared/{singularity_cache,lib}`.

## Tool versions per sample

Each per-sample process runs `bin/record_versions.sh` at the end of its script.
It runs the real `--version` commands inside the task container and writes
`<outdir>/sample_versions/<sample>/versions/<STEP>.versions.yml` (container image,
pipeline version/commit, host and time are included). Merge into one TSV per
sample, and optionally one combined TSV:

    bin/merge_versions.py results/sample_versions --all results/tool_versions.tsv

## AAFTF >= 0.7 command changes
This pipeline targets the AAFTF 0.7 command line. Changes that affect the modules:

- `-l/--left`, `-r/--right` are now `-1/--read1`, `-2/--read2` (`--merged` is unchanged).
- `--AAFTF_DB` was removed. AAFTF reads the database folder from the `AAFTF_DB`
  environment variable. The modules set `AAFTF_DB=/opt/aaftf_db`, the path the host
  database folder is bound to.
- Databases are downloaded only by `AAFTF database`. The host `AAFTF_DB` folder needs the
  new file `mitochondrion.1.1.genomic.fna.gz` for `vecscreen` (`AAFTF database mitodb`);
  the old `mito.nt.gz` is not used.
- `filter` on single-end reads writes `<prefix>_filtered_1` (it wrote `_filtered_U`). FILTER
  uses a separate prefix and renames the file.
- `fcs_gx_purge` has no `-c/--cpus` option, and it now drops only contigs that FCS-GX
  marks EXCLUDE (TRIM, FIX and REVIEW rows are kept with a warning).
- `assemble --method` is required, SPAdes runs with `--careful --isolate` by default, and
  `-m` is an integer number of GB.
- `AAFTF --version` / `-V` prints the git-tag version (for example `0.7.0-beta.4-15+g1f51deb`).

## Mitochondrial contig removal (MITOCHONDRION)

Runs after the contamination screens (FCS-GX, sourpurge) and before RMDUP. The earlier
screens do not look for organelle sequence: `fcs_screen`/`vecscreen` find adaptors and vectors
only, and FCS-GX reported no contigs for the MTK samples. Needs an AAFTF image with the depth MITO
flag (AAFTF main commit 0904010 or later). `aaftf_sif` is pinned to
`AAFTF-main-ccd6c26.sif` (AAFTF 0.7.0-beta.4-27+gccd6c26) until a stable 0.7.0 release.

**MITOCHONDRION** (on by default, `--skip_mitochondrion true` bypasses):
1. `AAFTF depth` maps the filtered reads to the assembly and flags candidate mitochondrial
   contigs: depth >= `mito_fold` (5) x the nuclear baseline, GC >= `mito_gc_drop` (10) points
   below the assembly, length >= `mito_min_len` (1000).
2. Only those candidates are searched against the RefSeq mitochondrial genomes. A candidate is
   removed only if the union of its alignments (>= `organelle_min_pident` 80% identity) covers
   >= `mito_min_cov` (0.5) of the contig. A NUMT inside a nuclear contig covers a small fraction
   and the contig is not a candidate anyway. `--mitochondrion_require_hit false` removes all
   candidates without this check.
3. rDNA veto: a contig that is >= `rdna_veto_cov` (0.2) rDNA (NCBI fungal and eukaryote-wide 18S/SSU, 28S/LSU, ITS) is never
   removed (`--rdna_veto_cov 0` disables the veto). On the MTK samples the 5 mitochondrial contigs
   have 0% rDNA coverage and 5 of the 6 other high-depth contigs (the rDNA repeat unit and its
   fragments) have 42-92% per database.
4. Removed contigs: `results/mitochondrion/<sample>.mitochondrion.fasta.gz`. Evidence:
   `<sample>.mitochondrion.{candidates,validation}.tsv` and `<sample>.mitochondrion.depth.txt`.
5. Kept but reported in `<sample>.mitochondrion.warnings.txt` and as a pipeline warning:
   `MITO_KEPT_UNVALIDATED` and `MITO_KEPT_RDNA` (candidates not removed) and `OUTLIER` (depth > mean +
   3 SD, not MITO), each OUTLIER labeled `rDNA-like` or `unexplained`. OUTLIER contigs are never removed.

**Plastid (chloroplast) contigs** are not handled here. A `PLASTID` step exists on the branch
`feat/plastid-step`. It is not in this release: on a cyanobacterial negative control it removed 4
rRNA-operon contigs (587-1,953 bp) that map at 100% identity to the genome, because cyanobacterial rRNA
matches cyanelle and chromatophore plastid genomes at 88-91% identity. It will come back after it is
tested on plastid-bearing data (a plant, a green alga, a lichen) with a depth gate.

**Databases.** Build once on a login node (compute nodes may have no internet):

    module load ncbi-blast/2.16.0+
    bin/download_organelle_db.sh all        # mitochondrion | rdna | all; -o DIR

This downloads the RefSeq mitochondrion `*.genomic.fna.gz` files (names are read from the directory
index) and NCBI's fungal and eukaryote-wide rRNA/ITS BLAST databases into `params.organelle_db`, and
builds the BLAST DBs there. Only accessions matching `ORGANELLE_ACC_REGEX` (default `^NC_`) are indexed
for mitochondrion. The RefSeq mitochondrion file also holds shotgun scaffolds (`NW_`, `NT_`),
and in release 237 `NW_026622727.1` (8.6 kb) matched non-mitochondrial 7 kb high-copy contigs at 94.5%
identity; with `^NC_` only, those contigs covered 9-27% and the mitochondrial contigs covered 100%.
The release number is stored in `RELEASE_NUMBER` and recorded in the version files.

## First-time setup

Two kinds of reference data are needed. The pipeline checks for them when it starts and stops within
seconds, printing the exact command to run, if something it needs is missing. It never downloads
anything by itself.

1. **AAFTF's contamination database (`contamdb`)**, about 5 MB. The pinned AAFTF (main `ccd6c26`) no longer
   reads the old `contamdb` text file; `FILTER` needs `contamdb.fa` and its BLAST files. Build it once per
   project, in the launch directory:

       nextflow run stajichlab/nf_aaftf -profile aaftf --setup_db true

   It writes to `params.aaftf_db`, default `<launchDir>/db`, and skips what already exists. Pass
   `--aaftf_db <folder>` to use another folder (for example one shared by several projects) for both
   `SETUP_DB` and the run. It downloads UniVec_Core, PhiX and the BBTools lambda and adapter sequences, and
   leaves any other file in the folder alone.
2. **Mitochondrion and rDNA databases** (about 700 MB, shared by all projects), built by
   `bin/download_organelle_db.sh` as described under *Mitochondrial contig removal* below.

`params.aaftf_db_shared` (default `null`) can name a second, read-only folder for large databases that many
projects share. Both folders are bound into the container and listed in `$AAFTF_DB` (first the project
folder, which AAFTF writes to, then the shared one).

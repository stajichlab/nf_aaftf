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
 ├─ ASSEMBLE ────────  spades (default) │ megahit │ unicycler │ dipspades
 │
 ├─ VECSCREEN ───────  AAFTF vecscreen  (optional, ON by default)
 │
 ├─ FCS_GX ──────────  NCBI FCS-GX purge (OPTIONAL, OFF by default)
 │
 ├─ RMDUP ───────────  remove duplicate contigs
 │
 ├─ POLISH ──────────  polca (default) │ pilon │ racon │ nextpolish
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
| ASSEMBLE | `assembler` | **spades** \| megahit \| unicycler \| dipspades | Only spades ships in the current image |
| ASSEMBLE | `assembler_args` | (string, default null) | Extra CLI args passed to the assembler |
| POLISH | `polisher` | **polca** \| pilon \| racon \| nextpolish | POLCA uses original filtered reads |

Outputs land under `results/<step>/` (publishDir), one file set per sample.

## Container

The pipeline is container-only. Set a single engine flag:

| `params.container_engine` | Container used |
|---|---|
| `singularity` (default) | cached SIF `AAFTF.sif` (from `ghcr.io/stajichlab/aaftf:latest`) |
| `docker`                 | `docker://ghcr.io/stajichlab/aaftf:latest` |

The AAFTF image bakes `AAFTF_DB=/opt/aaftf_db` and its toolchain
(`spades.py`, `bbduk.sh`, `fastp`, `pilon`, `polca.sh`, `minimap2`, `bwa`,
`bowtie2`, `taxonkit`, `sourmash` ...). The pipeline binds the host reference
DB and the NCBI taxdump into that path automatically via
`singularity.runOptions`.

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
| `assembler` | **spades** \| megahit \| unicycler \| dipspades | Assembly method |
| `assembler_args` | string (default null) | Extra CLI args for the assembler |
| `polisher` | **polca** \| pilon \| racon \| nextpolish | Polishing tool |

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
| `skip_vecscreen` | false | VecScreen on assembled contigs |
| `run_depth` | true | Read depth back-mapping |
| `fcs_taxid` | 4751 | Fallback NCBI taxonomy id for FCS-GX (4751 = Fungi) |

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

### FCS-GX (optional)
When enabled, `FCS_GX` resolves the **phylum** taxid at runtime from the
sample's organism taxid via `taxonkit` (mirroring the BFD reference
`GENOME_CLEAN` module) and stages the FCS-GX database into node-local scratch
before running `AAFTF fcs_gx_purge`. This step needs a current gxdb path
(`params.fcsgx_db`) and high memory; it is deliberately off by default.

### Alternative assemblers
- `assembler` is wired through to `AAFTF assemble --method`. Supported values:
  `spades` (default), `megahit`, `unicycler`, `dipspades`.
- **Only SPAdes ships in the current AAFTF image.** Megahit and unicycler require
  adding them to the container (or a future image build). The pipeline framework
  is ready; the tool availability is the only gap.
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
  (FCS_GX / SPAdes escalation). Verify with `sinfo`.
- AAFTF image + reference DBs are lab-shared under
  `/bigdata/stajichlab/shared/{singularity_cache,lib}`.

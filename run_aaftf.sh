#!/usr/bin/env bash
# run_aaftf.sh — sbatch launcher for the nf_aaftf Nextflow head process.
#
# Submit with:  sbatch run_aaftf.sh   [extra nextflow args...]
#
# The head process must outlive the whole run and must NOT live on the login
# node; it submits each AAFTF step as its own SLURM job through the slurm
# executor. Keep the head on a long-limit partition with a small footprint.

#SBATCH -p epyc
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 7-00:00:00
#SBATCH --job-name nf_aaftf
#SBATCH -o logs/slurm/nf_aaftf_%j.out
#SBATCH -e logs/slurm/nf_aaftf_%j.err

set -euo pipefail

mkdir -p logs/slurm logs/nextflow

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow
module load singularity

# Resolve the pipeline by PROJECT NAME, not this script's path (sbatch copies
# the script to a spool dir, so $0/$BASH_SOURCE are wrong).
#   PIPELINE  source of the pipeline (default: this project dir).
#             For a published project: PIPELINE=stajichlab/nf_aaftf
#   REVISION  git branch/tag/commit to run.
PIPELINE="${PIPELINE:-$PWD}"
REVISION="${REVISION:-}"

nextflow run "${PIPELINE}" ${REVISION:+-r "${REVISION}"} \
    -profile aaftf \
    -resume \
    "$@"

# Notes:
#  - Outputs and the sample sheet are read from the LAUNCH dir (${launchDir}),
#    i.e. where you sbatch from. Pass --samples /path/to/samples.csv if different.
#  - Cheap validation first (runs in seconds on a login node):
#        nextflow run "${PIPELINE}" -profile test -stub-run --n_test 2
#  - Choose the container engine with:  --container_engine docker   (or singularity)

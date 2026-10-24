// One-time setup of AAFTF's own contamination database (contamdb.fa + BLAST files).
// Run through the SETUP_DB entry workflow:
//   nextflow run <pipeline> -profile aaftf --setup_db true [--aaftf_db <folder>]
// It writes straight into params.aaftf_db (bound at /opt/aaftf_db), skips work that is already done,
// and is never cached by Nextflow, so deleting the folder and running it again rebuilds it.
process BUILD_CONTAMDB {
    tag   { 'contamdb' }
    label 'aaftf_lite'
    cache false

    script:
    """
    set -euo pipefail
    trap 'echo "SETUP_DB: ERROR: command failed at line \$LINENO of .command.sh (exit \$?)" >&2' ERR
    export AAFTF_DB=/opt/aaftf_db
    AAFTF database contamdb
    for f in contamdb.fa contamdb.nhr contamdb.nin contamdb.nsq; do
        [ -s "/opt/aaftf_db/\$f" ] || { echo "SETUP_DB: ERROR: \$f was not built (folder on the host: ${params.aaftf_db})" >&2; exit 1; }
    done
    echo "SETUP_DB: contamdb ready in ${params.aaftf_db}"
    """

    stub:
    """
    echo "SETUP_DB (stub): contamdb would be built in ${params.aaftf_db}"
    """
}

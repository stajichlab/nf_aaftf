#!/usr/bin/env bash
# download_organelle_db.sh — download RefSeq organelle genomes and build BLAST databases for the
# nf_aaftf ORGANELLE (mitochondrion) and PLASTID steps.
#
#   bin/download_organelle_db.sh [-o OUTDIR] [mitochondrion|rdna|all]
#     all = mitochondrion + rdna.  (A plastid database is not built: the PLASTID step is on feat/plastid-step.)
#
# OUTDIR defaults to /bigdata/stajichlab/shared/lib/AAFTF_DB/organelle (params.organelle_db).
# Creates OUTDIR/<kind>/{<kind>.*.genomic.fna.gz, RELEASE_NUMBER, FILES.txt, <kind>.n*} where the
# BLAST DB prefix is OUTDIR/<kind>/<kind>. The file list is read from the NCBI directory index
# (RefSeq may split a kind into several genomic files), so new RefSeq splits are picked up.
# rdna: NCBI's fungal 18S/28S/ITS and eukaryote-wide SSU/LSU/ITS BLAST databases
# (ftp.ncbi.nlm.nih.gov/blast/db), joined into one alias database OUTDIR/rdna/rdna. The pipeline uses it
# to protect nuclear rDNA contigs from removal (rDNA is high copy and can look like an organelle in the
# depth data). The eukaryote-wide sets cover plants and algae, which the fungal sets do not.
# Only curated RefSeq genomes (accessions matching ORGANELLE_ACC_REGEX, default ^NC_) go into the
# BLAST DB. The mitochondrion release also holds shotgun scaffolds (NW_, NT_), some of which are
# mislabeled or contain non-organelle sequence (release 237: NW_026622727.1, an 8.6 kb "mitochondrion"
# record that a nuclear 7 kb high-copy contig matched at 94.5% identity). Set ORGANELLE_ACC_REGEX='.'
# to keep everything. Only *.genomic.fna.gz is used; the .gbff and .protein.faa files are not needed for whole-contig
# screening. Run on a login node (compute nodes may lack internet) with makeblastdb on PATH:
#   module load ncbi-blast/2.16.0+ && bin/download_organelle_db.sh all
set -euo pipefail

BASE="https://ftp.ncbi.nlm.nih.gov/refseq/release"
OUT="/bigdata/stajichlab/shared/lib/AAFTF_DB/organelle"
while getopts "o:h" opt; do
    case "$opt" in
        o) OUT="$OPTARG" ;;
        h|*) sed -n '2,15p' "$0"; exit 0 ;;
    esac
done
shift $((OPTIND - 1))
WHAT="${1:-all}"
case "$WHAT" in mitochondrion|rdna) KINDS=("$WHAT") ;; all) KINDS=(mitochondrion rdna) ;; *) echo "unknown: $WHAT" >&2; exit 1 ;; esac
ACC_REGEX="${ORGANELLE_ACC_REGEX:-^NC_}"
command -v makeblastdb >/dev/null || { echo "ERROR: makeblastdb not on PATH (module load ncbi-blast/2.16.0+)" >&2; exit 1; }

release=$(curl -fsS --retry 3 "$BASE/RELEASE_NUMBER" | tr -d '[:space:]')
echo "RefSeq release: $release"
for kind in "${KINDS[@]}"; do
    dir="$OUT/$kind"; mkdir -p "$dir"
    if [ "$kind" = rdna ]; then
        command -v blastdb_aliastool >/dev/null || { echo "ERROR: blastdb_aliastool not on PATH" >&2; exit 1; }
        names="18S_fungal_sequences 28S_fungal_sequences ITS_RefSeq_Fungi SSU_eukaryote_rRNA LSU_eukaryote_rRNA ITS_eukaryote_sequences"
        for n in $names; do
            curl -fsS --retry 5 -o "$dir/$n.tar.gz" "https://ftp.ncbi.nlm.nih.gov/blast/db/$n.tar.gz"
            curl -fsS --retry 5 -o "$dir/$n.tar.gz.md5" "https://ftp.ncbi.nlm.nih.gov/blast/db/$n.tar.gz.md5"
            (cd "$dir" && md5sum -c "$n.tar.gz.md5" >/dev/null) || { echo "ERROR: md5 mismatch for $n" >&2; exit 1; }
            tar -xzf "$dir/$n.tar.gz" -C "$dir" && rm -f "$dir/$n.tar.gz"
        done
        [ -f "$dir/rdna.nal" ] && [ ! -f "$dir/rdna_fungal.nal" ] && cp -p "$dir/rdna.nal" "$dir/rdna_fungal.nal"   # keep the first (fungal-only) alias for rollback
        rm -f "$dir/rdna.nal"
        (cd "$dir" && blastdb_aliastool -dblist "$names" -dbtype nucl -out rdna -title "NCBI fungal + eukaryote 18S/SSU, 28S/LSU, ITS" >/dev/null)
        echo "$release" > "$dir/RELEASE_NUMBER"
        echo "rdna: alias DB at $dir/rdna ($names)"
        continue
    fi
    files=$(curl -fsS --retry 3 "$BASE/$kind/" | grep -oE "$kind\.[0-9]+\.[0-9]+\.genomic\.fna\.gz" | sort -u)
    [ -n "$files" ] || { echo "ERROR: no $kind genomic.fna.gz files found at $BASE/$kind/" >&2; exit 1; }
    echo "$kind: $(echo "$files" | wc -l) file(s): $(echo $files)"
    for f in $files; do
        # -C - resumes a partial download; the gzip test below catches a truncated file.
        curl -fsS --retry 5 -C - -o "$dir/$f" "$BASE/$kind/$f"
        gzip -t "$dir/$f" || { echo "ERROR: $dir/$f failed gzip -t" >&2; exit 1; }
    done
    echo "$files" > "$dir/FILES.txt"
    # Build into a temp prefix then move, so a failed build never leaves a half DB in place.
    tmp="$dir/.build_$$"; mkdir -p "$tmp"
    (cd "$dir" && zcat $files) \
        | awk -v re="$ACC_REGEX" '/^>/{keep = (substr($1,2) ~ re)} keep' \
        | makeblastdb -in - -dbtype nucl -parse_seqids -title "RefSeq $kind release $release (accessions $ACC_REGEX)" -out "$tmp/$kind" >/dev/null
    rm -f "$dir/$kind".n*
    mv "$tmp/$kind".n* "$dir/"; rmdir "$tmp"
    echo "$release" > "$dir/RELEASE_NUMBER"
    echo "$kind: BLAST DB at $dir/$kind (release $release, $(blastdbcmd -db "$dir/$kind" -info | grep -E 'sequences' | head -1 | xargs))"
done

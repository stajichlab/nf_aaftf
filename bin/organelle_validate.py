#!/usr/bin/env python3
"""Validate candidate organelle contigs by whole-contig coverage against an organelle BLAST DB.

Only the candidate contigs are searched (never the whole assembly), and a contig is validated
only when the union of its alignments covers at least --min-cov of its length at >= --min-pident
identity. A nuclear contig carrying a short NUMT/NUPT insertion covers a small fraction and fails.

usage: organelle_validate.py --candidates FASTA --db BLASTDB --out TSV --validated IDS
           [--min-cov 0.5] [--min-pident 80] [--evalue 1e-10] [--cpus 4]
Writes TSV columns: contig, length, covered_bp, covered_frac, best_pident, n_hsps, top_subject, validated
and IDS (one validated contig id per line). Needs blastn on PATH; standard library only.
"""
import argparse
import subprocess
import sys


def read_lengths(fasta):
    lengths, name = {}, None
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name:
                lengths[name] += len(line.strip())
    return lengths


def union_bp(intervals):
    total, end = 0, 0
    for s, e in sorted(intervals):
        if e <= end:
            continue
        total += e - max(s, end)
        end = e
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--validated", required=True)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--min-pident", type=float, default=80.0)
    ap.add_argument("--evalue", default="1e-10")
    ap.add_argument("--cpus", type=int, default=4)
    a = ap.parse_args()

    lengths = read_lengths(a.candidates)
    hits = {}
    if lengths:
        cmd = ["blastn", "-task", "blastn", "-query", a.candidates, "-db", a.db, "-evalue", a.evalue,
               "-num_threads", str(a.cpus), "-max_target_seqs", "50",
               "-outfmt", "6 qseqid qstart qend pident sseqid bitscore"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            sys.exit("ERROR: blastn failed: " + res.stderr.strip())
        for line in res.stdout.splitlines():
            q, qs, qe, pid, sid, bits = line.split("\t")
            pid = float(pid)
            if pid < a.min_pident:
                continue
            qs, qe = sorted((int(qs), int(qe)))
            hits.setdefault(q, []).append((qs - 1, qe, pid, sid, float(bits)))

    validated = []
    with open(a.out, "w") as out:
        out.write("contig\tlength\tcovered_bp\tcovered_frac\tbest_pident\tn_hsps\ttop_subject\tvalidated\n")
        for c, length in lengths.items():
            h = hits.get(c, [])
            cov = union_bp([(s, e) for s, e, *_ in h])
            frac = cov / length if length else 0.0
            ok = frac >= a.min_cov
            best = max((x[2] for x in h), default=0.0)
            top = max(h, key=lambda x: x[4])[3] if h else "-"
            out.write(f"{c}\t{length}\t{cov}\t{frac:.3f}\t{best:.1f}\t{len(h)}\t{top}\t{'yes' if ok else 'no'}\n")
            if ok:
                validated.append(c)
    with open(a.validated, "w") as fh:
        fh.writelines(c + "\n" for c in validated)


if __name__ == "__main__":
    main()

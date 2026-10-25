#!/usr/bin/env python3
"""Decide which organelle (mitochondrion or plastid) each candidate contig belongs to, if any.

A candidate is scored by the fraction of its length covered by alignments to RefSeq mitochondrial genomes
(--mito), plastid genomes (--plast) and rDNA (--rdna); each table is organelle_validate.py output. A table that
is not given counts as zero coverage, so without the plastid database a candidate can only be called
mitochondrial, which is how the pipeline behaved before plastids were added.

Order of the rules for one contig:
 1. rDNA coverage >= --rdna-veto-cov (and the veto is on): kept, never removed.
 2. The database with the higher coverage wins. If the two differ by less than --tie-eps, the contig is called
    mitochondrial only when the depth/GC (MITO) rule flagged it (AT-rich, high-copy); otherwise plastid.
    Plant mitochondrial genomes carry chloroplast DNA, so a chloroplast contig can cover 100% of both.
 3. The winner must cover at least its own threshold (--mito-min-cov, --plastid-min-cov), else the contig is
    kept as unvalidated.
 4. A contig called mitochondrial or plastid is removed only if that removal is switched on
    (--remove-mito, --remove-plastid); otherwise it is kept and reported.

usage: organelle_classify.py --candidates TSV --pool IDS --out TSV --remove-mito-ids FILE --remove-plastid-ids FILE
           --warnings FILE [--mito TSV] [--plast TSV] [--rdna TSV] [--assume-mito] [--remove-mito] [--remove-plastid]
           [--mito-min-cov 0.7] [--plastid-min-cov 0.4] [--rdna-veto-cov 0.2] [--tie-eps 0.1]
Standard library only.
"""
import argparse
import csv


def read_cov(path):
    """contig -> (covered_frac, best_pident, top_subject) from an organelle_validate.py table; {} if no path."""
    out = {}
    if not path:
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["contig"]] = (float(r["covered_frac"]), float(r["best_pident"]), r["top_subject"])
    return out


def classify(rule, cm, cp, cr, mito_min_cov, plastid_min_cov, rdna_veto_cov, tie_eps, remove_mito, remove_plastid):
    """Return (class, decision) for one contig. class is mito, plastid or '-'."""
    if rdna_veto_cov > 0 and cr >= rdna_veto_cov:
        return "-", "kept_rdna"
    if max(cm, cp) <= 0:
        return "-", "kept_unvalidated"
    if cm - cp > tie_eps:
        cls = "mito"
    elif cp - cm > tie_eps:
        cls = "plastid"
    else:
        cls = "mito" if rule else "plastid"
    cov, need = (cm, mito_min_cov) if cls == "mito" else (cp, plastid_min_cov)
    if cov < need:
        return cls, "kept_unvalidated"
    if cls == "mito":
        return cls, "removed_mito" if remove_mito else "kept_mito_candidate"
    return cls, "removed_plastid" if remove_plastid else "kept_plastid_candidate"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", required=True, help="AAFTF depth mito_candidates.tsv (every contig)")
    ap.add_argument("--pool", required=True, help="file of contig ids to classify")
    ap.add_argument("--out", required=True)
    ap.add_argument("--remove-mito-ids", required=True)
    ap.add_argument("--remove-plastid-ids", required=True)
    ap.add_argument("--warnings", required=True)
    ap.add_argument("--mito")
    ap.add_argument("--plast")
    ap.add_argument("--rdna")
    ap.add_argument("--assume-mito", action="store_true", help="no mito table: treat MITO-rule candidates as fully covered")
    ap.add_argument("--remove-mito", action="store_true")
    ap.add_argument("--remove-plastid", action="store_true")
    ap.add_argument("--mito-min-cov", type=float, default=0.7)
    ap.add_argument("--plastid-min-cov", type=float, default=0.4)
    ap.add_argument("--rdna-veto-cov", type=float, default=0.2)
    ap.add_argument("--tie-eps", type=float, default=0.1)
    a = ap.parse_args()

    info = {}
    with open(a.candidates) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            info[r["contig"]] = r
    pool = [l.strip() for l in open(a.pool) if l.strip()]
    mito, plast, rdna = read_cov(a.mito), read_cov(a.plast), read_cov(a.rdna)

    rm_mito, rm_plast, warn = [], [], []
    with open(a.out, "w") as out:
        out.write("contig\tlength\tfold\tgc_percent\tmito_rule\tmito_cov\tmito_pident\tmito_top\t"
                  "plast_cov\tplast_pident\tplast_top\trdna_frac\tclass\tdecision\n")
        for c in pool:
            r = info[c]
            rule = r["mito"] == "yes"
            cm, pm, tm = mito.get(c, (0.0, 0.0, "-"))
            if a.assume_mito and not a.mito and rule:
                cm, pm, tm = 1.0, 0.0, "assumed"
            cp, pp, tp = plast.get(c, (0.0, 0.0, "-"))
            cr = rdna.get(c, (0.0, 0.0, "-"))[0]
            cls, dec = classify(rule, cm, cp, cr, a.mito_min_cov, a.plastid_min_cov, a.rdna_veto_cov, a.tie_eps,
                                a.remove_mito, a.remove_plastid)
            out.write("\t".join(map(str, [c, r["length"], r["fold_over_baseline"], r["gc_percent"], "yes" if rule else "no",
                                          f"{cm:.3f}", f"{pm:.1f}", tm, f"{cp:.3f}", f"{pp:.1f}", tp, f"{cr:.3f}", cls, dec])) + "\n")
            if dec == "removed_mito":
                rm_mito.append(c)
            elif dec == "removed_plastid":
                rm_plast.append(c)
            else:
                origin = "MITO" if rule else "PLASTID"
                warn.append(f"{origin}_{dec.upper()}\t{c}\tlen={r['length']}\tmito_cov={cm:.2f}\tplast_cov={cp:.2f}\trdna_cov={cr:.2f}")
    for path, ids in ((a.remove_mito_ids, rm_mito), (a.remove_plastid_ids, rm_plast)):
        with open(path, "w") as fh:
            fh.writelines(i + "\n" for i in ids)
    with open(a.warnings, "w") as fh:
        fh.writelines(w + "\n" for w in warn)


if __name__ == "__main__":
    main()

"""Unit tests for bin/organelle_classify.py. Run: python3 -m unittest tests/test_organelle_classify.py -v"""
import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "bin", "organelle_classify.py")
spec = importlib.util.spec_from_file_location("organelle_classify", SCRIPT)
oc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oc)


def run(rule, cm, cp, cr=0.0, rm=True, rp=True, mito_min=0.7, plast_min=0.4, veto=0.2, eps=0.1):
    return oc.classify(rule, cm, cp, cr, mito_min, plast_min, veto, eps, rm, rp)


class Classify(unittest.TestCase):
    def test_rdna_veto_beats_everything(self):
        self.assertEqual(run(True, 1.0, 0.0, cr=0.82), ("-", "kept_rdna"))
        self.assertEqual(run(False, 0.0, 1.0, cr=0.2), ("-", "kept_rdna"))

    def test_veto_off(self):
        self.assertEqual(run(True, 1.0, 0.0, cr=0.9, veto=0), ("mito", "removed_mito"))

    def test_clear_mito(self):
        self.assertEqual(run(True, 0.94, 0.02), ("mito", "removed_mito"))

    def test_clear_plastid(self):
        self.assertEqual(run(False, 0.01, 0.99), ("plastid", "removed_plastid"))

    def test_chloroplast_contig_that_matches_a_plant_mitogenome(self):
        # Chlamydomonas NODE_4712: 53% mito (lavender mitogenome with chloroplast DNA), 99% plastid, AT-rich so the MITO rule fired
        self.assertEqual(run(True, 0.53, 0.99), ("plastid", "removed_plastid"))

    def test_tie_goes_to_mito_only_if_the_mito_rule_fired(self):
        # Arabidopsis chloroplast contigs: 100% in both databases
        self.assertEqual(run(False, 1.0, 1.0), ("plastid", "removed_plastid"))
        self.assertEqual(run(True, 1.0, 1.0), ("mito", "removed_mito"))

    def test_thresholds_are_per_class(self):
        self.assertEqual(run(True, 0.69, 0.0), ("mito", "kept_unvalidated"))      # relatives-only mitogenome, just under 0.7
        self.assertEqual(run(False, 0.0, 0.39), ("plastid", "kept_unvalidated"))
        self.assertEqual(run(False, 0.0, 0.40), ("plastid", "removed_plastid"))
        self.assertEqual(run(True, 0.5, 0.0, mito_min=0.5), ("mito", "removed_mito"))

    def test_without_the_plastid_database_only_mito_is_possible(self):
        self.assertEqual(run(True, 0.8, 0.0), ("mito", "removed_mito"))
        self.assertEqual(run(False, 0.8, 0.0), ("mito", "removed_mito"))   # coverage alone decides; callers only pass candidates

    def test_no_hit(self):
        self.assertEqual(run(True, 0.0, 0.0), ("-", "kept_unvalidated"))

    def test_removal_switched_off_keeps_and_reports(self):
        self.assertEqual(run(True, 0.9, 0.0, rm=False), ("mito", "kept_mito_candidate"))
        self.assertEqual(run(False, 0.0, 0.9, rp=False), ("plastid", "kept_plastid_candidate"))


class Cli(unittest.TestCase):
    def test_end_to_end(self):
        with tempfile.TemporaryDirectory() as d:
            p = lambda n: os.path.join(d, n)
            open(p("cand.tsv"), "w").write(
                "contig\tlength\tmean_depth\tfold_over_baseline\tgc_percent\tmito\n"
                "m1\t15000\t300\t30\t33\tyes\n" "c1\t5000\t700\t60\t45\tyes\n" "c2\t8000\t900\t80\t37\tno\n"
                "r1\t2000\t500\t40\t50\tno\n" "n1\t3000\t12\t1\t55\tno\n")
            hdr = "contig\tlength\tcovered_bp\tcovered_frac\tbest_pident\tn_hsps\ttop_subject\tvalidated\n"
            open(p("mito.tsv"), "w").write(hdr + "m1\t15000\t14000\t0.93\t99\t3\tNC_1\tyes\n" "c1\t5000\t2650\t0.53\t100\t2\tNC_2\tyes\n"
                                           "c2\t8000\t8000\t1.0\t100\t1\tNC_3\tyes\n")
            open(p("plast.tsv"), "w").write(hdr + "c1\t5000\t4950\t0.99\t98\t2\tNC_9\tyes\n" "c2\t8000\t8000\t1.0\t100\t1\tNC_8\tyes\n")
            open(p("rdna.tsv"), "w").write(hdr + "r1\t2000\t1700\t0.85\t90\t1\tX\tno\n")
            open(p("pool.ids"), "w").write("m1\nc1\nc2\nr1\n")
            cmd = [sys.executable, SCRIPT, "--candidates", p("cand.tsv"), "--pool", p("pool.ids"), "--mito", p("mito.tsv"),
                   "--plast", p("plast.tsv"), "--rdna", p("rdna.tsv"), "--out", p("out.tsv"), "--remove-mito-ids", p("rm.ids"),
                   "--remove-plastid-ids", p("rp.ids"), "--warnings", p("warn.txt"), "--remove-mito", "--remove-plastid"]
            subprocess.run(cmd, check=True)
            self.assertEqual(open(p("rm.ids")).read().split(), ["m1"])
            # c1: 53% mito vs 99% plastid -> plastid; c2: ties at 1.0 and the MITO rule did not fire -> plastid
            self.assertEqual(open(p("rp.ids")).read().split(), ["c1", "c2"])
            self.assertIn("PLASTID_KEPT_RDNA\tr1", open(p("warn.txt")).read())


if __name__ == "__main__":
    unittest.main()

import csv
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from benchmarks.rvne_evaluation import average_precision, boundary_metrics, comparator_profiles, correlations
from data_generation.rvne_simulation import SCENARIOS, simulate
from msa_conservation.cli import main
from msa_conservation.msa import Alignment
from msa_conservation.scoring.rvne import RVNEParameters, score_reference, site_profiles


class IntegrationTests(unittest.TestCase):
    def test_cli_writes_scores_and_segments(self):
        with TemporaryDirectory() as root:
            folder = Path(root)
            msa = folder / "tiny.fasta"
            msa.write_text(">a\nAARRRVA\n>b\nAARRRLA\n>c\nRAR--VA\n", encoding="utf-8")
            sites, segments = folder / "sites.tsv", folder / "segments.tsv"
            result = main(["score", "--input", str(msa), "--method", "scaw-rvne",
                           "--output", str(sites), "--segments-output", str(segments)])
            self.assertEqual(result, 0)
            with sites.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 7)
            self.assertEqual(rows[3]["position"], "4")
            self.assertEqual(rows[3]["score"], "")
            self.assertEqual(rows[0]["r"], "0.8")
            with segments.open(newline="", encoding="utf-8") as handle:
                regions = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(regions[0]["start"], "1")
            self.assertEqual(regions[-1]["end"], "7")

    def test_synthetic_families_are_reproducible(self):
        for name in SCENARIOS:
            first = simulate(scenario=name, size=20, length=100, seed=70)
            second = simulate(scenario=name, size=20, length=100, seed=70)
            self.assertEqual(first.alignment, second.alignment)
            self.assertEqual(first.truth.shape, (100,))
            np.testing.assert_array_equal(first.truth, second.truth)
        motif = simulate(scenario="motifs", size=20, motif_length=3)
        self.assertEqual(motif.conserved.sum(), 9)
        self.assertEqual(len(motif.boundaries), 6)

    def test_supplemental_patterns_have_the_requested_profiles(self):
        homogeneous = simulate(scenario="conserved", length=24)
        blocks = simulate(scenario="conserved-variable-blocks", length=24, variable_fraction=0.25)
        motif = simulate(scenario="single-motif", length=24, motif_length=3)

        np.testing.assert_array_equal(homogeneous.rates, [0.06] * 24)
        self.assertEqual(homogeneous.boundaries, ())
        np.testing.assert_array_equal(blocks.rates, [0.06] * 9 + [1.0] * 6 + [0.06] * 9)
        self.assertEqual(blocks.boundaries, (9, 15))
        np.testing.assert_array_equal(motif.rates, [1.0] * 10 + [0.06] * 3 + [1.0] * 11)
        self.assertEqual(motif.boundaries, (10, 13))

    def test_three_rate_pattern_recovers_both_flanks_at_odd_length(self):
        pattern = simulate(scenario="three-rate-blocks", length=25, divergence=0.7)
        expected_rates = np.array([0.06] * 4 + [0.32] * 3 + [1.0] * 10 + [0.32] * 3 + [0.06] * 5)

        np.testing.assert_array_equal(pattern.rates, expected_rates)
        np.testing.assert_allclose(pattern.truth, np.exp(-0.7 * expected_rates))
        np.testing.assert_array_equal(pattern.conserved, expected_rates == 0.06)
        self.assertEqual(pattern.boundaries, (4, 7, 17, 20))

    def test_evaluation_ties_and_boundary_matching(self):
        scores = np.array([3., 3., 2., 1.])
        self.assertAlmostEqual(correlations(scores, scores, rank=True), 1.)
        self.assertAlmostEqual(average_precision(scores, np.array([True, False, True, False])),
                               (.5 + 2/3) / 2)
        self.assertEqual(boundary_metrics((10, 11, 22), (11, 20)), (2/3, 1., .8))

    def test_comparator_consensus_coverage_and_kl_equivalence(self):
        alignment = Alignment(("a", "b", "c", "d"),
                              ("AAA-AXAA", "ASV-RXA-", "ATAAAZAA", "ATAAA-AA"))
        sites = site_profiles(alignment)
        scores = comparator_profiles(sites, alignment)
        np.testing.assert_allclose(scores["ClustalW consensus"], [3, 2, 1, np.nan, 0, np.nan, 3, 0])
        expected_kl = score_reference(sites, RVNEParameters(1, 0, 0)) * sites.coverage
        np.testing.assert_allclose(scores["KL"], expected_kl * np.log2(20), atol=1e-12)

    def test_al2co_variance_uses_pooled_unweighted_frequencies(self):
        alignment = Alignment(("a", "b", "c", "d"), ("AS", "AS", "AS", "T-"))
        scores = comparator_profiles(site_profiles(alignment), alignment)
        pooled = np.array([3, 3, 1]) / 7
        expected = np.square(np.array([[.75, 0, .25], [0, 1, 0]]) - pooled).sum(axis=1)
        np.testing.assert_allclose(scores["AL2CO variance"], expected * [1, .75])

import unittest
from io import StringIO
from math import fsum, log

from msa_conservation.msa import AA_INDEX, Alignment, henikoff_weights, read_fasta
from msa_conservation.scoring.fusion import score_alignment
from msa_conservation.scoring.substitution import load_blosum62


class BaselineScoringTests(unittest.TestCase):
    def test_missing_columns_and_coverage(self):
        alignment = read_fasta(StringIO(">one description\na-XA\n>two\na.X-\n"))
        conserved, gap, unknown, partial = tuple(score_alignment(alignment, gap_penalty=0.4))
        self.assertEqual(alignment.ids, ("one", "two"))
        self.assertEqual(conserved.entropy_variability, 0.0)
        self.assertEqual(conserved.substitution_variability, 0.0)
        self.assertEqual(conserved.integrated_variability, 0.0)
        self.assertEqual(gap.nongap_coverage, 0.0)
        self.assertEqual(unknown.nongap_coverage, 1.0)
        self.assertEqual(unknown.canonical_coverage, 0.0)
        self.assertIsNone(gap.entropy_variability)
        self.assertIsNone(unknown.integrated_variability)
        self.assertAlmostEqual(partial.integrated_variability, 0.2)

    def test_blosum_cost_and_fusion_use_the_same_column_frequencies(self):
        alignment = Alignment(("one", "two"), ("AI", "WL"))
        radical, conservative = tuple(score_alignment(alignment))
        matrix = load_blosum62()
        cost = matrix.severity[AA_INDEX["A"]][AA_INDEX["W"]] / matrix.normalizer
        self.assertAlmostEqual(fsum(matrix.background), 1.0)
        self.assertAlmostEqual(radical.entropy_variability, log(2) / log(20))
        self.assertEqual(radical.entropy_variability, conservative.entropy_variability)
        self.assertAlmostEqual(radical.conditional_severity, cost)
        self.assertAlmostEqual(radical.substitution_variability, cost / 2)
        self.assertAlmostEqual(radical.integrated_variability, log(2) / log(20) * (1 + cost))
        self.assertGreater(radical.integrated_variability, conservative.integrated_variability)
        no_fusion = tuple(score_alignment(alignment, gamma=0))
        self.assertEqual(no_fusion[0].entropy_variability, no_fusion[0].integrated_variability)

    def test_henikoff_downweights_redundancy(self):
        alignment = Alignment(("a", "b", "c"), ("A", "A", "R"))
        weights = henikoff_weights(alignment)
        self.assertEqual(weights, (0.25, 0.25, 0.5))
        site = next(score_alignment(alignment, weights=weights))
        self.assertAlmostEqual(site.entropy_variability, log(2) / log(20))

    def test_reject_unaligned_input_and_zero_weights(self):
        with self.assertRaisesRegex(ValueError, "expected 2"):
            read_fasta(StringIO(">one\nAA\n>two\nA\n"))
        alignment = Alignment(("a",), ("A",))
        with self.assertRaisesRegex(ValueError, "positive"):
            next(score_alignment(alignment, weights=(0.0,)))

import unittest
from io import StringIO
from math import log, log2

from msa_conservation.msa import AMINO_ACIDS, Alignment, henikoff_weights, read_fasta
from msa_conservation.scoring.entropy import adjusted_entropy, jensen_shannon_divergence, relative_entropy
from msa_conservation.scoring.baseline import score_alignment


class BaselineScoringTests(unittest.TestCase):
    def test_missing_columns_and_coverage(self):
        alignment = read_fasta(StringIO(">one description\na-XA\n>two\na.X-\n"))
        conserved, gap, unknown, partial = tuple(score_alignment(alignment, gap_penalty=0.4))
        self.assertEqual(alignment.ids, ("one", "two"))
        self.assertEqual(conserved.entropy_variability, 0.0)
        self.assertEqual(gap.nongap_coverage, 0.0)
        self.assertEqual(unknown.nongap_coverage, 1.0)
        self.assertEqual(unknown.canonical_coverage, 0.0)
        self.assertIsNone(gap.entropy_variability)
        self.assertIsNone(unknown.entropy_variability)
        self.assertAlmostEqual(partial.entropy_variability, 0.2)

    def test_shannon_depends_only_on_residue_frequencies(self):
        alignment = Alignment(("one", "two"), ("AI", "WL"))
        first, second = score_alignment(alignment)
        self.assertAlmostEqual(first.entropy_variability, log(2) / log(20))
        self.assertEqual(first.entropy_variability, second.entropy_variability)

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

    def test_relative_entropy_and_jsd(self):
        background = (0.25, 0.75)
        self.assertEqual(relative_entropy(background, background), 0)
        self.assertEqual(relative_entropy((1, 0), background), log2(4))
        self.assertEqual(jensen_shannon_divergence(background, background), 0)
        self.assertEqual(jensen_shannon_divergence((1, 0), (0, 1)), 1)

    def test_adjusted_entropy_uses_unweighted_distinct_residue_pairs(self):
        # The published examples use rounded row means. Use their exact integer sums.
        for residues, observed_mean, background_mean in (("PFN", -3, -94 / 57), ("YWF", 2, -97 / 57)):
            frequencies = tuple(1 / 3 if residue in residues else 0 for residue in AMINO_ACIDS)
            coefficient = 2 ** ((observed_mean - background_mean) / 2)
            self.assertAlmostEqual(adjusted_entropy(frequencies), log(3) / coefficient)
        self.assertAlmostEqual(adjusted_entropy((1 / 20,) * 20), log(20))

    def test_weighted_sum_of_pairs_excludes_self_comparisons_and_gaps(self):
        alignment = Alignment(("a1", "a2", "r", "gap"), ("AA-", "AA-", "R--", "-XX"))
        first, second, third = score_alignment(alignment, weights=(1, 2, 3, 4), reference_baselines=True)
        # Distinct pairs at the first site contribute 1*2*4 + 1*3*(-1) + 2*3*(-1).
        self.assertAlmostEqual(first.blosum62_sop_variability, 1 / 11)
        self.assertAlmostEqual(second.blosum62_sop_variability, -4)
        self.assertEqual(second.adjusted_entropy_variability, 0)
        self.assertIsNone(third.blosum62_sop_variability)
        self.assertIsNone(third.jsd_variability)
        singleton = next(score_alignment(Alignment(("a", "b"), ("A", "-")), reference_baselines=True))
        self.assertIsNone(singleton.blosum62_sop_variability)
        self.assertLess(singleton.jsd_variability, 0)

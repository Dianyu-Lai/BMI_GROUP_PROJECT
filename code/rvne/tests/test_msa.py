import unittest
from io import StringIO

from msa_conservation.msa import AA_INDEX, Alignment, column_profiles, henikoff_weights, read_fasta


class AlignmentTests(unittest.TestCase):
    def test_missing_columns_and_coverage(self):
        alignment = read_fasta(StringIO(">one description\na-XA\n>two\na.X-\n"))
        conserved, gap, unknown, partial = column_profiles(alignment)
        self.assertEqual(alignment.ids, ("one", "two"))
        self.assertEqual(conserved.frequencies[AA_INDEX["A"]], 1.0)
        self.assertEqual(gap.nongap_coverage, 0.0)
        self.assertEqual(unknown.nongap_coverage, 1.0)
        self.assertEqual(unknown.canonical_coverage, 0.0)
        self.assertEqual(partial.canonical_coverage, 0.5)

    def test_henikoff_downweights_redundancy(self):
        alignment = Alignment(("a", "b", "c"), ("A", "A", "R"))
        weights = henikoff_weights(alignment)
        self.assertEqual(weights, (0.25, 0.25, 0.5))
        site = next(column_profiles(alignment, weights=weights))
        self.assertEqual(site.frequencies[AA_INDEX["A"]], 0.5)
        self.assertEqual(site.frequencies[AA_INDEX["R"]], 0.5)

    def test_reject_unaligned_input_and_zero_weights(self):
        with self.assertRaisesRegex(ValueError, "expected 2"):
            read_fasta(StringIO(">one\nAA\n>two\nA\n"))
        alignment = Alignment(("a",), ("A",))
        with self.assertRaisesRegex(ValueError, "positive"):
            next(column_profiles(alignment, weights=(0.0,)))

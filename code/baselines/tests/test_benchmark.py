import unittest

from benchmarks.metrics import (
    matched_pair_accuracy,
    roc_auc,
    select_threshold,
    sensitivity_specificity,
    spearman,
)
from benchmarks.synthetic import RESIDUE_CLASSES, generate_alignment, generate_constraint_alignment
from msa_conservation.scoring.fusion import score_alignment


class BenchmarkTests(unittest.TestCase):
    def test_progressive_divergence_preserves_ancestry_and_truth(self):
        parameters = dict(size=16, length=90, seed=4)
        ancestor = generate_alignment(**parameters, divergence=0.0)
        shallow = generate_alignment(**parameters, divergence=0.3)
        deep = generate_alignment(**parameters, divergence=1.0)
        self.assertEqual(ancestor.site_rates, deep.site_rates)
        self.assertEqual(sum(deep.variable_sites), 30)
        for root, low, high in zip(
            ancestor.alignment.sequences, shallow.alignment.sequences, deep.alignment.sequences
        ):
            for root_residue, low_residue, high_residue in zip(root, low, high):
                if root_residue != low_residue:
                    self.assertEqual(low_residue, high_residue)

    def test_tie_aware_metrics_and_threshold(self):
        scores = [0.0, 0.5, 0.5, 1.0]
        labels = [False, False, True, True]
        self.assertEqual(roc_auc(scores, labels), 0.875)
        threshold = select_threshold(scores, labels)
        self.assertEqual(threshold, 1.0)
        self.assertEqual(sensitivity_specificity(scores, labels, threshold), (0.5, 1.0))
        self.assertEqual(spearman(scores, scores), 1.0)
        self.assertIsNone(spearman([1.0] * len(scores), scores))

    def test_constraint_pairs_match_mutations_and_entropy(self):
        parameters = dict(size=16, length=90, seed=4)
        roots = generate_constraint_alignment(**parameters, divergence=0.0)
        dataset = generate_constraint_alignment(**parameters, divergence=0.8)
        self.assertEqual(roots.variable_sites, dataset.variable_sites)
        self.assertEqual(len(set(dataset.site_rates)), 1)
        classes = {residue: group for group in RESIDUE_CLASSES for residue in group}
        for root, sequence in zip(roots.alignment.sequences, dataset.alignment.sequences):
            for i in range(0, dataset.alignment.length, 2):
                self.assertEqual(sequence[i] == root[i], sequence[i + 1] == root[i + 1])
            for ancestral, residue, label in zip(root, sequence, dataset.variable_sites):
                if ancestral != residue:
                    self.assertEqual(classes[ancestral] != classes[residue], label)

        entropy = [site.entropy_variability for site in score_alignment(dataset.alignment)]
        self.assertEqual(roc_auc(entropy, dataset.variable_sites), 0.5)
        self.assertEqual(matched_pair_accuracy(entropy, dataset.variable_sites), 0.5)
        self.assertEqual(matched_pair_accuracy((0.9, 0.1), (True, False)), 1.0)
        self.assertEqual(matched_pair_accuracy((0.9, 0.1), (False, True)), 0.0)

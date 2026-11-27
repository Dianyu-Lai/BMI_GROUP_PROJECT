import unittest
from math import exp, fsum, log
from statistics import mean, pstdev

from msa_conservation.scoring.physicochemical import (
    aaindex_features,
    aaindex_squared_distances,
    calibrated_blosum62,
)
from msa_conservation.scoring.substitution import load_blosum62


class PhysicochemicalTests(unittest.TestCase):
    def test_standardization_and_eta_zero(self):
        for values in zip(*aaindex_features()):
            self.assertAlmostEqual(mean(values), 0.0)
            self.assertAlmostEqual(pstdev(values), 1.0)
        self.assertIs(calibrated_blosum62(0), load_blosum62())

    def test_tilt_and_diagonal_normalized_correction(self):
        baseline = load_blosum62()
        distances = aaindex_squared_distances()
        eta = 0.5
        matrix = calibrated_blosum62(eta)
        z = fsum(
            p * exp(-eta * distance / 2)
            for row, dists in zip(baseline.joint, distances)
            for p, distance in zip(row, dists)
        )
        self.assertAlmostEqual(fsum(matrix.background), 1.0)
        self.assertEqual(matrix.normalizer, baseline.normalizer)
        for i, row in enumerate(matrix.joint):
            self.assertAlmostEqual(matrix.background[i], fsum(row))
            for j, probability in enumerate(row):
                expected_q = baseline.joint[i][j] * exp(-eta * distances[i][j] / 2) / z
                self.assertAlmostEqual(probability, expected_q)
                self.assertEqual(matrix.severity[i][j], matrix.severity[j][i])
                log_odds = baseline.log_odds
                base_cost = (log_odds[i][i] + log_odds[j][j]) / 2 - log_odds[i][j]
                expected_cost = max(0, base_cost + eta * distances[i][j] / log(2))
                self.assertAlmostEqual(matrix.severity[i][j], expected_cost)

    def test_invalid_tilt(self):
        for eta in (-1, float("inf"), float("nan")):
            with self.subTest(eta=eta), self.assertRaisesRegex(ValueError, "finite and nonnegative"):
                calibrated_blosum62(eta)

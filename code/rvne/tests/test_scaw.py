import itertools
import unittest

import numpy as np

from msa_conservation.scoring.scaw import (
    ContextParameters, WeightedIntervals, context_from_settings, contextual_scores,
    noise_variance, segmentation_penalty, segment_profile,
)


class ScawTests(unittest.TestCase):
    def test_pelt_matches_exhaustive_search(self):
        rng = np.random.default_rng(812)
        for size in range(1, 9):
            for multiplier in (0., .1, 1., 4.):
                for _ in range(6):
                    values = rng.normal(size=size)
                    weights = rng.uniform(.1, 1., size=size)
                    intervals = WeightedIntervals(values, weights)
                    penalty = segmentation_penalty(intervals, ContextParameters(multiplier))
                    segments = segment_profile(intervals, ContextParameters(multiplier))
                    actual = sum(intervals.moments(s.start, s.stop)[2] for s in segments)
                    actual += penalty * (len(segments) - 1)
                    expected = np.inf
                    for n in range(size):
                        for boundaries in itertools.combinations(range(1, size), n):
                            points = (0, *boundaries, size)
                            cost = sum(intervals.moments(a, b)[2]
                                       for a, b in zip(points, points[1:])) + penalty * n
                            expected = min(expected, cost)
                    self.assertAlmostEqual(actual, expected, places=10)

    def test_windows_do_not_cross_short_motif(self):
        raw = np.r_[np.zeros(20), np.ones(3), np.zeros(20)]
        coverage = np.ones(len(raw))
        result = contextual_scores(raw, coverage, "scaw")
        self.assertEqual([s.stop for s in result.segments], [20, 23, 43])
        np.testing.assert_allclose(result.scores, raw, atol=1e-14)
        self.assertEqual(result.window_size[21], 3)
        self.assertEqual(result.window_size[10], 7)
        self.assertLess(contextual_scores(raw, coverage, "uniform").scores[21], .5)

    def test_missing_and_tiny_profiles(self):
        raw = np.array([.3, np.nan, .3, np.nan])
        coverage = np.array([1., .4, .9, 0.])
        for method in ("uniform", "gaussian", "segmentation", "scaw"):
            result = contextual_scores(raw, coverage, method)
            np.testing.assert_allclose(result.scores, [.3, np.nan, .27, np.nan])
            tiny = contextual_scores(raw[:1], coverage[:1], method)
            np.testing.assert_allclose(tiny.scores, [.3])
            missing = contextual_scores(raw[1:2], coverage[1:2], method)
            self.assertTrue(np.isnan(missing.scores[0]))
            empty = contextual_scores(raw[:0], coverage[:0], method)
            self.assertEqual(len(empty.scores), 0)

    def test_noise_estimator_and_penalty_units(self):
        raw, weights = np.array([1., 2., 1., 3.]), np.array([1., .5, .8, 1.])
        intervals = WeightedIntervals(raw, weights)
        noise = np.mean(np.diff(raw)**2 / (1 / weights[:-1] + 1 / weights[1:]))
        self.assertAlmostEqual(noise_variance(intervals), noise)
        self.assertAlmostEqual(segmentation_penalty(intervals, ContextParameters()),
                               noise * 2 * np.log(4))
        unit = WeightedIntervals(raw, np.ones(4))
        self.assertAlmostEqual(noise_variance(unit), sum(np.diff(raw)**2) / 6)

    def test_missing_runs_are_not_compressed_for_noise(self):
        raw = np.array([np.nan, .1, np.nan, .9, np.nan])
        weights = np.array([0., 1., .2, 1., 0.])
        intervals = WeightedIntervals(raw, weights)
        self.assertEqual(noise_variance(intervals), intervals.variance)
        parts = segment_profile(intervals, ContextParameters(0))
        self.assertEqual([(s.start, s.stop) for s in parts], [(0, 3), (3, 5)])
        result = contextual_scores(raw, weights, "scaw", ContextParameters(0))
        np.testing.assert_allclose(result.scores, raw)

    def test_constant_ties_and_affine_invariance(self):
        rng = np.random.default_rng(77)
        raw = np.r_[rng.normal(0, .1, 30), rng.normal(1, .1, 30)]
        weights = rng.uniform(.7, 1., len(raw))
        first = contextual_scores(raw, weights, "scaw")
        second = contextual_scores(3 * raw + 2, weights, "scaw")
        self.assertEqual([(s.start, s.stop) for s in first.segments],
                         [(s.start, s.stop) for s in second.segments])
        np.testing.assert_allclose(second.estimates, 3 * first.estimates + 2)
        for value in (0., .3, 1e6):
            parts = segment_profile(WeightedIntervals(np.full(20, value), weights[:20]),
                                    ContextParameters(0))
            self.assertEqual(len(parts), 1)

    def test_fixed_center_neighbor_formula(self):
        raw = np.array([0., 2., 4.])
        coverage = np.array([1., .5, 1.])
        result = contextual_scores(raw, coverage, "scaw", ContextParameters(100))
        np.testing.assert_allclose(result.estimates, [(.5 * (1 + 4) / 1.5), 2.,
                                                     2 + .5 * 1 / 1.5])
        np.testing.assert_allclose(result.scores, result.estimates * coverage)

    def test_invalid_inputs_and_obsolete_settings(self):
        for value in (-1, float("nan"), float("inf"), True, "1"):
            with self.assertRaises(ValueError):
                ContextParameters(value)
        with self.assertRaisesRegex(ValueError, "v1 settings"):
            context_from_settings({"beta": 6, "smoothing": .25})
        for values, weights in (([1], [2]), ([1], [float("nan")]),
                                ([float("inf")], [1]), ([1, 2], [1])):
            with self.assertRaises(ValueError):
                WeightedIntervals(np.array(values), np.array(weights))

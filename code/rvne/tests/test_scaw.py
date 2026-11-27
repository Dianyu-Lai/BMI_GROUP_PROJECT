import itertools
import unittest

import numpy as np

from msa_conservation.scoring.scaw import (
    ContextParameters, WeightedIntervals, contextual_scores, segment_profile,
)


class ScawTests(unittest.TestCase):
    def test_segment_dp_matches_exhaustive_search(self):
        values = np.array([.1, .2, .1, .8, .8, .9, .3, .3])
        weights = np.array([1., .7, .9, 1., .8, .8, .9, 1.])
        params = ContextParameters(beta=.4, min_segment=2, max_segment=5)
        intervals = WeightedIntervals(values, weights)
        penalty = params.beta * intervals.variance
        segments = segment_profile(intervals, params)
        actual = sum(intervals.moments(s.start, s.stop)[2] + penalty for s in segments)
        expected = np.inf
        for n in range(len(values)):
            for boundaries in itertools.combinations(range(1, len(values)), n):
                points = (0, *boundaries, len(values))
                if all(2 <= b-a <= 5 for a, b in zip(points, points[1:])):
                    expected = min(expected, sum(intervals.moments(a, b)[2] + penalty
                                                 for a, b in zip(points, points[1:])))
        self.assertAlmostEqual(actual, expected)

    def test_windows_do_not_cross_short_motif(self):
        raw = np.r_[np.zeros(20), np.ones(3), np.zeros(20)]
        coverage = np.ones(len(raw))
        result = contextual_scores(raw, coverage, "scaw", ContextParameters(beta=2))
        self.assertEqual([s.stop for s in result.segments], [20, 23, 43])
        np.testing.assert_allclose(result.scores, raw, atol=1e-14)
        self.assertEqual(result.window_size[21], 3)
        self.assertGreater(result.window_size[10], result.window_size[21])
        self.assertLess(contextual_scores(raw, coverage, "uniform").scores[21], .3)

    def test_missing_and_tiny_profiles(self):
        raw = np.array([.3, np.nan, .3, np.nan])
        coverage = np.array([1., .4, .9, 0.])
        for method in ("uniform", "gaussian", "adaptive", "segmentation", "scaw"):
            result = contextual_scores(raw, coverage, method)
            np.testing.assert_allclose(result.scores, [.3, np.nan, .27, np.nan])
            tiny = contextual_scores(raw[:1], coverage[:1], method)
            np.testing.assert_allclose(tiny.scores, [.3])
            missing = contextual_scores(raw[1:2], coverage[1:2], method)
            self.assertTrue(np.isnan(missing.scores[0]))

import unittest

import numpy as np

from msa_conservation.msa import Alignment
from msa_conservation.scoring.rvne import RVNEParameters, score_reference, site_profiles
from msa_conservation.scoring.rvne_optimized import score_grid, score_optimized


class OptimizedRVNETests(unittest.TestCase):
    def test_reference_grid_and_optimized_agree(self):
        alignment = Alignment(("one", "two", "three", "four"),
                              ("ARL-DW", "ASQ-DW", "TRC-DM", "TRC-DX"))
        sites = site_profiles(alignment)
        grid = [RVNEParameters(r, s1, s2)
                for r in (0, .8, 1) for s1 in (0, .8, 1) for s2 in (0, .8, .99)]
        result = score_grid(sites, grid, batch_size=2)
        for i, params in enumerate(grid):
            expected = score_reference(sites, params)
            np.testing.assert_allclose(score_optimized(sites, params, batch_size=2),
                                       expected, atol=1e-12, rtol=0)
            np.testing.assert_allclose(result[i], expected[sites.valid], atol=1e-12, rtol=0)

    def test_all_missing(self):
        sites = site_profiles(Alignment(("one",), ("-X",)))
        params = RVNEParameters()
        self.assertTrue(np.isnan(score_optimized(sites, params)).all())
        self.assertEqual(score_grid(sites, [params]).shape, (1, 0))

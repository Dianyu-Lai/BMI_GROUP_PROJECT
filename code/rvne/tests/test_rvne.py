import unittest

import numpy as np

from msa_conservation.msa import AA_INDEX, Alignment, henikoff_weights
from msa_conservation.scoring.rvne import (
    DEFAULT_PARAMETERS, SIZE, RVNEParameters, background, blosum_pair_frequencies,
    density, paper_parameter_grid, relative_entropy, relative_entropy_trace,
    score_reference, similarity, site_profiles,
)


class RVNETests(unittest.TestCase):
    def test_reference_frequencies_and_paper_table(self):
        m = blosum_pair_frequencies()
        self.assertAlmostEqual(m.sum(), 0.9987)
        np.testing.assert_array_equal(m, m.T)
        indices = [AA_INDEX[aa] for aa in "LMW"]
        conditional = m[np.ix_(indices, indices)]
        np.testing.assert_allclose(conditional / conditional.sum(),
                                   [[.627, .083, .012], [.083, .068, .003],
                                    [.012, .003, .110]], atol=.0005)

    def test_similarity_endpoints_and_psd(self):
        np.testing.assert_array_equal(similarity(0), np.eye(SIZE))
        np.testing.assert_array_equal(similarity(1), np.ones((SIZE, SIZE)))
        for s in [i / 10 for i in range(1, 10)] + [.99]:
            matrix = similarity(s)
            np.testing.assert_allclose(matrix.diagonal(), 1)
            self.assertGreater(np.linalg.eigvalsh(matrix)[0], 0)

    def test_background_and_grid(self):
        uniform = np.full(SIZE, 1 / SIZE)
        np.testing.assert_allclose(background(0), uniform)
        self.assertAlmostEqual(background(1).sum(), 1)
        np.testing.assert_allclose(background(.8), .2 * uniform + .8 * background(1))
        grid = paper_parameter_grid()
        self.assertEqual(len(set(grid)), 1331)
        self.assertIn(DEFAULT_PARAMETERS, grid)
        with self.assertRaises(ValueError):
            RVNEParameters(s2=1)

    def test_identical_and_redundant_sequences(self):
        alignment = Alignment(("a", "b", "c"), ("AAA", "AAA", "RRR"))
        np.testing.assert_allclose(henikoff_weights(alignment), (.25, .25, .5))
        identical = Alignment(("a", "b"), ("ARN", "ARN"))
        np.testing.assert_allclose(henikoff_weights(identical), (.5, .5))
        sites = site_profiles(alignment)
        np.testing.assert_allclose(sites.probabilities[:, :2], .5)

    def test_gap_cutoff_and_ambiguity(self):
        # First column exactly 30% gaps, second 40%, third has no canonical observations.
        seqs = ("AAX",) * 6 + ("A-X", "-.X", "-.X", "-.X")
        sites = site_profiles(Alignment(tuple(map(str, range(10))), seqs))
        np.testing.assert_allclose(sites.coverage, [.7, .6, 1])
        np.testing.assert_array_equal(sites.valid, [True, False, False])
        np.testing.assert_allclose(sites.probabilities[:2, 0], 1)
        raw = score_reference(sites)
        self.assertGreater(raw[0] * sites.coverage[0], 0)
        self.assertTrue(np.isnan(raw[1:]).all())

    def test_density_and_special_cases(self):
        p = np.array([.4, .6] + [0.] * (SIZE - 2))
        q = background(1)
        diagonal, mixed = density(p, similarity(0)), density(p, similarity(.8))
        for matrix in (diagonal, mixed, density(q, similarity(.99))):
            np.testing.assert_allclose(matrix, matrix.T)
            self.assertAlmostEqual(np.trace(matrix), 1)
            self.assertGreaterEqual(np.linalg.eigvalsh(matrix)[0], -1e-12)
        uniform = np.eye(SIZE) / SIZE
        p_positive = p[:2]
        self.assertAlmostEqual(relative_entropy(diagonal, uniform),
                               1 + (p_positive @ np.log(p_positive)) / np.log(SIZE))
        eigenvalues = np.linalg.eigvalsh(mixed)[-2:]
        self.assertAlmostEqual(relative_entropy(mixed, uniform),
                               1 + (eigenvalues @ np.log(eigenvalues)) / np.log(SIZE))
        self.assertAlmostEqual(relative_entropy(diagonal, np.diag(q)),
                               (p_positive @ np.log(p_positive / q[:2])) / np.log(SIZE))
        pure = density(p, similarity(1))
        self.assertAlmostEqual(relative_entropy(pure, uniform), 1)
        with self.assertRaises(ValueError):
            relative_entropy(mixed, np.diag(np.r_[np.zeros(19), 1.]))

    def test_reference_and_independent_trace_form(self):
        rng = np.random.default_rng(217)
        for s1 in (0, .8, 1):
            for s2 in (0, .8, .99):
                p = rng.dirichlet(np.ones(SIZE))
                rho = density(p, similarity(s1))
                sigma = density(background(.8), similarity(s2))
                self.assertAlmostEqual(relative_entropy(rho, sigma),
                                       relative_entropy_trace(rho, sigma), places=12)
                self.assertAlmostEqual(relative_entropy(sigma, sigma), 0, places=12)


if __name__ == "__main__":
    unittest.main()

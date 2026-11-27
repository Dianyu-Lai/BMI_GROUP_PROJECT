"""Cached, batched trace formulation and reuse across RVNE parameter settings."""

from functools import lru_cache
from typing import Sequence

import numpy as np

from .rvne import (
    DEFAULT_PARAMETERS, RVNEParameters, SiteProfiles, _background_eigenpairs,
    _spectral_entropy_terms, background, density, similarity,
)


@lru_cache(maxsize=192)
def log_background(r: float, s2: float) -> np.ndarray:
    """Cache the base-20 spectral logarithm of the fixed background density."""
    log_values, vectors = _background_eigenpairs(density(background(r), similarity(s2)))
    result = (vectors * log_values) @ vectors.T
    result.setflags(write=False)
    return result


def score_optimized(
    sites: SiteProfiles, parameters: RVNEParameters = DEFAULT_PARAMETERS, *, batch_size: int = 256,
) -> np.ndarray:
    """Return reference-equivalent raw RVNE using observed eigenvalues only."""
    if batch_size < 1:
        raise ValueError("batch size must be positive")
    selected = np.flatnonzero(sites.valid)
    result = np.full(len(sites.valid), np.nan)
    observed = similarity(parameters.s1)
    log_sigma = log_background(parameters.r, parameters.s2)
    for start in range(0, len(selected), batch_size):
        indices = selected[start:start + batch_size]
        rho = density(sites.probabilities[indices], observed)
        terms = _spectral_entropy_terms(np.linalg.eigvalsh(rho))
        result[indices] = terms - np.einsum("ijk,jk->i", rho, log_sigma)
    return result


def score_grid(
    sites: SiteProfiles, grid: Sequence[RVNEParameters], *, batch_size: int = 256,
) -> np.ndarray:
    """Score all parameters on usable columns, in grid order and then site order."""
    if batch_size < 1:
        raise ValueError("batch size must be positive")
    probabilities = sites.probabilities[sites.valid]
    result = np.empty((len(grid), len(probabilities)))
    for s1 in sorted({p.s1 for p in grid}):
        position = [i for i, p in enumerate(grid) if p.s1 == s1]
        observed = similarity(s1)
        log_sigmas = np.stack([log_background(grid[i].r, grid[i].s2) for i in position])
        kernels = (log_sigmas * observed).reshape(len(position), -1)
        for start in range(0, len(probabilities), batch_size):
            stop = start + batch_size
            roots = np.sqrt(probabilities[start:stop])
            outer = roots[:, :, None] * roots[:, None, :]
            rho = outer * observed
            terms = _spectral_entropy_terms(np.linalg.eigvalsh(rho))
            result[position, start:stop] = terms - kernels @ outer.reshape(len(roots), -1).T
    return result

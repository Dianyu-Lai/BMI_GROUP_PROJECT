"""Johansson and Toh's relative von Neumann entropy, with an independent trace check."""

from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from itertools import product
from math import log

import numpy as np

from ..msa import Alignment, AMINO_ACIDS, GAPS, column_profiles, henikoff_weights


SIZE = len(AMINO_ACIDS)
LOG_SIZE = log(SIZE)
PSD_TOLERANCE = 1e-12
TENTHS = tuple(i / 10 for i in range(11))


@dataclass(frozen=True, order=True)
class RVNEParameters:
    """Background relativity, observed similarity and background similarity."""

    r: float = 0.8
    s1: float = 0.8
    s2: float = 0.8

    def __post_init__(self) -> None:
        if not (0 <= self.r <= 1 and 0 <= self.s1 <= 1 and 0 <= self.s2 < 1):
            raise ValueError("require 0 <= r,s1 <= 1 and 0 <= s2 < 1")


DEFAULT_PARAMETERS = RVNEParameters()


@dataclass(frozen=True)
class SiteProfiles:
    """Column probabilities and unweighted nongap fraction, in alignment order."""

    probabilities: np.ndarray
    coverage: np.ndarray
    valid: np.ndarray


@lru_cache(maxsize=1)
def blosum_pair_frequencies() -> np.ndarray:
    """Read the symmetric joint frequencies from the rounded original qij table."""
    source = files("msa_conservation").joinpath("data", "blosum62.qij")
    lines = [line.split() for line in source.read_text(encoding="ascii").splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    if "".join(lines[0]) != AMINO_ACIDS or len(lines[1:]) != SIZE:
        raise ValueError("unexpected BLOSUM62 qij alphabet or matrix size")
    matrix = np.zeros((SIZE, SIZE))
    for i, row in enumerate(lines[1:]):
        if len(row) != i + 1:
            raise ValueError("expected a lower-triangular qij table")
        matrix[i, :i + 1] = matrix[:i + 1, i] = tuple(map(float, row))
    if not np.all(np.isfinite(matrix) & (matrix > 0)) or abs(matrix.sum() - 1) > 0.02:
        raise ValueError("invalid BLOSUM62 qij probabilities")
    matrix.setflags(write=False)
    return matrix


@lru_cache(maxsize=64)
def similarity(s: float) -> np.ndarray:
    """Construct A(s), verifying positive semidefiniteness with a symmetric solver."""
    if not 0 <= s <= 1:
        raise ValueError("similarity must be in [0, 1]")
    if s == 0:
        matrix = np.eye(SIZE)
    elif s == 1:
        matrix = np.ones((SIZE, SIZE))
    else:
        joint = blosum_pair_frequencies()
        normalized = joint / np.sqrt(np.outer(joint.diagonal(), joint.diagonal()))
        matrix = normalized ** ((1 - s) / s)
    if np.linalg.eigvalsh(matrix)[0] < -PSD_TOLERANCE:
        raise ValueError("A(s) is not positive semidefinite")
    matrix.setflags(write=False)
    return matrix


def background(r: float) -> np.ndarray:
    """Interpolate between the uniform and normalized BLOSUM62 marginals."""
    if not 0 <= r <= 1:
        raise ValueError("background relativity must be in [0, 1]")
    joint = blosum_pair_frequencies()
    return (1 - r) / SIZE + r * joint.sum(axis=1) / joint.sum()


def density(probabilities: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    """Construct one density matrix or a batch from normalized probabilities."""
    roots = np.sqrt(probabilities)
    return roots[..., :, None] * matrix * roots[..., None, :]


def _spectral_entropy_terms(values: np.ndarray) -> np.ndarray:
    """Evaluate sum λ log_20(λ) along the last axis."""
    if np.any(values < -PSD_TOLERANCE):
        raise ValueError("observed density has a materially negative eigenvalue")
    positive = values > 0
    logs = np.log(values, out=np.zeros_like(values), where=positive)
    return np.sum(values * logs, axis=-1) / LOG_SIZE


def _background_eigenpairs(sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eigenvalues, vectors = np.linalg.eigh(sigma)
    if np.any(eigenvalues <= 0):
        raise ValueError("background density must be numerically positive definite")
    return np.log(eigenvalues) / LOG_SIZE, vectors


def relative_entropy(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Calculate the original eigenvector-overlap expression (equation 11)."""
    values, vectors = np.linalg.eigh(rho)
    log_background, background_vectors = _background_eigenpairs(sigma)
    overlaps = (vectors.T @ background_vectors) ** 2
    return float(_spectral_entropy_terms(values) - np.maximum(values, 0) @ (overlaps @ log_background))


def relative_entropy_trace(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Calculate the equivalent matrix trace expression for baseline validation."""
    log_values, vectors = _background_eigenpairs(sigma)
    log_sigma = (vectors * log_values) @ vectors.T
    return float(_spectral_entropy_terms(np.linalg.eigvalsh(rho)) - np.trace(rho @ log_sigma))


def site_profiles(alignment: Alignment, max_gap: float = 0.3) -> SiteProfiles:
    """Calculate Henikoff-weighted frequencies and the paper's gap mask."""
    if not 0 <= max_gap <= 1:
        raise ValueError("maximum gap fraction must be in [0, 1]")
    weights = henikoff_weights(alignment)
    probabilities = np.array([site.frequencies for site in column_profiles(alignment, weights)])
    gap_fractions = np.fromiter(
        (sum(residue in GAPS for residue in column) / alignment.size
         for column in zip(*alignment.sequences)),
        dtype=float, count=alignment.length,
    )
    valid = (gap_fractions <= max_gap) & (probabilities.sum(axis=1) > 0)
    return SiteProfiles(probabilities, 1 - gap_fractions, valid)


def score_reference(
    sites: SiteProfiles, parameters: RVNEParameters = DEFAULT_PARAMETERS
) -> np.ndarray:
    """Return raw base-20 RVNE, marking excluded columns NaN."""
    observed = similarity(parameters.s1)
    sigma = density(background(parameters.r), similarity(parameters.s2))
    scores = np.full(len(sites.valid), np.nan)
    for i in np.flatnonzero(sites.valid):
        scores[i] = relative_entropy(density(sites.probabilities[i], observed), sigma)
    return scores


def paper_parameter_grid() -> tuple[RVNEParameters, ...]:
    """Return the 1,331 settings from the published search."""
    return tuple(RVNEParameters(*values) for values in product(TENTHS, TENTHS, (*TENTHS[:-1], 0.99)))

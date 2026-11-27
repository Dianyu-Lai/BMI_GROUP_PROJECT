"""Entropy and divergence measures for canonical amino acid frequencies."""

from functools import lru_cache
from itertools import combinations
from math import fsum, inf, log, log2
from statistics import mean
from typing import Sequence

from ..msa import AMINO_ACIDS
from .substitution import blosum62_integer_scores


LOG_ALPHABET_SIZE = log(len(AMINO_ACIDS))


def shannon_entropy(frequencies: Sequence[float]) -> float:
    """Return Shannon variability on a fixed 20-residue alphabet."""
    value = fsum(-p * log(p) for p in frequencies if p > 0) / LOG_ALPHABET_SIZE
    return min(1.0, value)


def relative_entropy(frequencies: Sequence[float], background: Sequence[float]) -> float:
    """Return the KL divergence in bits."""
    terms = (
        p * (log2(p) - log2(q)) if q > 0 else inf
        for p, q in zip(frequencies, background, strict=True)
        if p > 0
    )
    return max(0.0, fsum(terms))


def jensen_shannon_divergence(frequencies: Sequence[float], background: Sequence[float]) -> float:
    """Return the equal-weight Jensen–Shannon divergence in bits."""
    midpoint = tuple((p + q) / 2 for p, q in zip(frequencies, background, strict=True))
    return (relative_entropy(frequencies, midpoint) + relative_entropy(background, midpoint)) / 2


@lru_cache(maxsize=1)
def _blosum62_mean_off_diagonal() -> tuple[float, ...]:
    return tuple(
        mean(score for j, score in enumerate(row) if i != j)
        for i, row in enumerate(blosum62_integer_scores())
    )


def adjusted_entropy(frequencies: Sequence[float]) -> float:
    """Return Zhang et al.'s Shannon entropy divided by a similarity coefficient."""
    observed = [i for i, frequency in enumerate(frequencies) if frequency > 0]
    if len(observed) < 2:
        return 0.0
    matrix = blosum62_integer_scores()
    row_means = _blosum62_mean_off_diagonal()
    observed_similarity = mean(matrix[i][j] for i, j in combinations(observed, 2))
    expected_similarity = mean(row_means[i] for i in observed)
    coefficient = 2 ** ((observed_similarity - expected_similarity) / 2)
    entropy = fsum(-frequencies[i] * log(frequencies[i]) for i in observed)
    return entropy / coefficient

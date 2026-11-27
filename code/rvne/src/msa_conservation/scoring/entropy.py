"""Normalized Shannon entropy for canonical amino acid frequencies."""

from math import fsum, log
from typing import Sequence

from ..msa import AMINO_ACIDS


LOG_ALPHABET_SIZE = log(len(AMINO_ACIDS))


def shannon_entropy(frequencies: Sequence[float]) -> float:
    """Return Shannon variability on a fixed 20-residue alphabet."""
    value = fsum(-p * log(p) for p in frequencies if p > 0) / LOG_ALPHABET_SIZE
    return min(1.0, value)

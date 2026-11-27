"""BLOSUM62-based substitution severity from full-precision pair probabilities."""

import csv
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import fsum, isclose, log2
from pathlib import Path
from typing import Sequence

from ..msa import AMINO_ACIDS


REFERENCE_PATH = Path(__file__).resolve().parents[3] / "data" / "reference" / "blosum62_joint.tsv"
SIZE = len(AMINO_ACIDS)


@dataclass(frozen=True)
class Blosum62:
    joint: tuple[tuple[float, ...], ...]
    background: tuple[float, ...]
    severity: tuple[tuple[float, ...], ...]
    normalizer: float
    log_odds: tuple[tuple[float, ...], ...]


@lru_cache(maxsize=1)
def load_blosum62() -> Blosum62:
    with REFERENCE_PATH.open(encoding="ascii", newline="") as file:
        header, *rows = csv.reader(file, delimiter="\t")
    alphabet = list(AMINO_ACIDS)
    if header != ["aa", *alphabet] or [row[0] for row in rows] != alphabet:
        raise ValueError("BLOSUM62 table uses an unexpected amino acid order")
    joint = tuple(tuple(map(float, row[1:])) for row in rows)
    if any(len(row) != SIZE or any(not 0 < p < 1 for p in row) for row in joint):
        raise ValueError("BLOSUM62 table must contain 20 by 20 positive probabilities")
    if not isclose(fsum(p for row in joint for p in row), 1.0, rel_tol=1e-12):
        raise ValueError("BLOSUM62 probabilities must sum to one")
    if any(joint[i][j] != joint[j][i] for i, j in combinations(range(SIZE), 2)):
        raise ValueError("BLOSUM62 pair probabilities must be symmetric")

    background = tuple(fsum(row) for row in joint)
    log_odds = tuple(
        tuple(2 * log2(joint[i][j] / (background[i] * background[j])) for j in range(SIZE))
        for i in range(SIZE)
    )
    return _build_blosum62(joint, background, log_odds, normalizer=None)


def blosum_from_log_joint(
    log_joint: tuple[tuple[float, ...], ...], *, normalizer: float
) -> Blosum62:
    """Construct a matrix from a normalized, symmetric log2 joint distribution."""
    joint = tuple(tuple(2 ** log_p for log_p in row) for row in log_joint)
    background = tuple(fsum(row) for row in joint)
    log_background = tuple(log2(p) for p in background)
    log_odds = tuple(
        tuple(2 * (log_joint[i][j] - (log_background[i] + log_background[j])) for j in range(SIZE))
        for i in range(SIZE)
    )
    return _build_blosum62(joint, background, log_odds, normalizer=normalizer)


def _build_blosum62(
    joint: tuple[tuple[float, ...], ...],
    background: tuple[float, ...],
    log_odds: tuple[tuple[float, ...], ...],
    *,
    normalizer: float | None,
) -> Blosum62:
    severity = tuple(
        tuple(max(0.0, (log_odds[i][i] + log_odds[j][j]) / 2 - log_odds[i][j]) for j in range(SIZE))
        for i in range(SIZE)
    )
    if normalizer is None:
        unlike_pairs = tuple(combinations(range(SIZE), 2))
        unlike_mass = fsum(joint[i][j] for i, j in unlike_pairs)
        normalizer = fsum(joint[i][j] * severity[i][j] for i, j in unlike_pairs) / unlike_mass
    return Blosum62(joint, background, severity, normalizer, log_odds)


def substitution_scores(
    frequencies: Sequence[float], matrix: Blosum62
) -> tuple[float, float]:
    """Return P (expected cost) and M (cost conditional on unlike draws)."""
    observed = [(i, frequency) for i, frequency in enumerate(frequencies) if frequency > 0]
    probability_mass = []
    weighted_costs = []
    for (i, p), (j, q) in combinations(observed, 2):
        pair_mass = 2 * p * q  # Count both ordered pairs.
        probability_mass.append(pair_mass)
        weighted_costs.append(pair_mass * matrix.severity[i][j])

    unlike_mass = fsum(probability_mass)
    if unlike_mass == 0:
        return 0.0, 0.0
    expected_cost = fsum(weighted_costs) / matrix.normalizer
    return expected_cost, expected_cost / unlike_mass


@lru_cache(maxsize=1)
def blosum62_integer_scores() -> tuple[tuple[int, ...], ...]:
    """Return the canonical BLOSUM62 scores in integer half-bit units."""
    return tuple(tuple(round(score) for score in row) for row in load_blosum62().log_odds)


def blosum62_sum_of_pairs(
    frequencies: Sequence[float], same_residue_pair_mass: Sequence[float]
) -> float | None:
    """Return mean BLOSUM62 similarity over distinct pairs of sequences."""
    matrix = blosum62_integer_scores()
    observed = [(i, frequency) for i, frequency in enumerate(frequencies) if frequency > 0]
    masses = [same_residue_pair_mass[i] for i, _ in observed]
    scores = [mass * matrix[i][i] for mass, (i, _) in zip(masses, observed)]
    for (i, p), (j, q) in combinations(observed, 2):
        pair_mass = p * q
        masses.append(pair_mass)
        scores.append(pair_mass * matrix[i][j])
    total_mass = fsum(masses)
    return fsum(scores) / total_mass if total_mass else None

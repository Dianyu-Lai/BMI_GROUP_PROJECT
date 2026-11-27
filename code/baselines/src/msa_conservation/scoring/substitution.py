"""BLOSUM62 background frequencies, integer scores, and weighted sum of pairs."""

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
    background: tuple[float, ...]
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
    return Blosum62(background, log_odds)


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

"""BLOSUM62-based substitution severity from full-precision pair probabilities."""

import csv
from functools import lru_cache
from itertools import combinations
from math import fsum, isclose, log2
from pathlib import Path

from ..msa import AMINO_ACIDS


REFERENCE_PATH = Path(__file__).resolve().parents[3] / "data" / "reference" / "blosum62_joint.tsv"
SIZE = len(AMINO_ACIDS)


@lru_cache(maxsize=1)
def blosum62_severity() -> tuple[tuple[float, ...], ...]:
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
    return tuple(
        tuple(max(0.0, (log_odds[i][i] + log_odds[j][j]) / 2 - log_odds[i][j]) for j in range(SIZE))
        for i in range(SIZE)
    )

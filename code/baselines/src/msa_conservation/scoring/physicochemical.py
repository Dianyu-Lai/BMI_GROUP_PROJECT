"""AAindex1 residue descriptors and an RBF tilt of BLOSUM62 pair probabilities."""

import csv
from functools import lru_cache
from math import fsum, isfinite, log, log2
from statistics import mean, pstdev

from ..msa import AMINO_ACIDS
from .substitution import REFERENCE_PATH, Blosum62, blosum_from_log_joint, load_blosum62


PROPERTY_IDS = ("KYTJ820101", "GRAR740103", "KLEP840101", "VINM940101")
PROPERTY_PATH = REFERENCE_PATH.with_name("aaindex1_properties.tsv")


@lru_cache(maxsize=1)
def aaindex_features() -> tuple[tuple[float, ...], ...]:
    """Return one vector per residue, standardized over the 20 canonical residues."""
    with PROPERTY_PATH.open(encoding="ascii", newline="") as file:
        header, *rows = csv.reader(file, delimiter="\t")
    if header != ["aa", *PROPERTY_IDS] or [row[0] for row in rows] != list(AMINO_ACIDS):
        raise ValueError("AAindex table uses an unexpected feature or residue order")
    values = [tuple(map(float, row[1:])) for row in rows]
    if any(len(row) != len(PROPERTY_IDS) or not all(map(isfinite, row)) for row in values):
        raise ValueError("AAindex descriptors must be complete and finite")

    standardized = []
    for property_values in zip(*values):
        center = mean(property_values)
        scale = pstdev(property_values)
        if scale == 0:
            raise ValueError("AAindex descriptors must have nonzero variance")
        standardized.append(tuple((value - center) / scale for value in property_values))
    return tuple(zip(*standardized))


@lru_cache(maxsize=1)
def aaindex_squared_distances() -> tuple[tuple[float, ...], ...]:
    features = aaindex_features()
    return tuple(
        tuple(
            fsum((a - b) ** 2 for a, b in zip(first, second)) / len(PROPERTY_IDS)
            for second in features
        )
        for first in features
    )


@lru_cache(maxsize=16)
def calibrated_blosum62(eta: float) -> Blosum62:
    """Apply q_eta(a,b) proportional to q_0(a,b) exp(-eta d²(a,b)/2)."""
    if not isfinite(eta) or eta < 0:
        raise ValueError("AAindex eta must be finite and nonnegative")
    baseline = load_blosum62()
    if eta == 0:
        return baseline

    distances = aaindex_squared_distances()
    log_kernel_scale = eta / (2 * log(2))
    log_unnormalized = tuple(
        tuple(log2(p) - log_kernel_scale * distance for p, distance in zip(row, dists))
        for row, dists in zip(baseline.joint, distances)
    )
    # The diagonal retains positive mass even when distant off-diagonal entries underflow.
    log_total = log2(fsum(2 ** value for row in log_unnormalized for value in row))
    log_joint = tuple(tuple(value - log_total for value in row) for row in log_unnormalized)
    if any(not isfinite(2 * value) for row in log_joint for value in row):
        raise ValueError("AAindex eta is too large for finite substitution scores")
    return blosum_from_log_joint(log_joint, normalizer=baseline.normalizer)

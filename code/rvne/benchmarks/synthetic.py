"""Generate aligned protein sequences with controlled mutation or substitution labels."""

from dataclasses import dataclass
from math import expm1, isfinite
from random import Random

from msa_conservation.msa import AMINO_ACIDS, Alignment


# Low, intermediate, and high mutation rates in consecutive alignment regions.
REGION_RATES = (0.05, 0.25, 1.0)
RESIDUE_CLASSES = ("AVLIM", "FWY", "STNQ", "KR", "DE")
MATCHED_SITE_RATE = 1.0


@dataclass(frozen=True)
class SyntheticAlignment:
    alignment: Alignment
    site_rates: tuple[float, ...]
    variable_sites: tuple[bool, ...]


def generate_alignment(
    *, size: int, length: int, divergence: float, seed: int
) -> SyntheticAlignment:
    """Simulate independent descendants of one reference with at most one change per site."""
    if size < 2 or length < 3:
        raise ValueError("synthetic alignments need at least two sequences and three sites")
    _check_divergence(divergence)

    rng = Random(seed)
    alphabet_size = len(AMINO_ACIDS)
    ancestor = [rng.randrange(alphabet_size) for _ in range(length)]
    site_rates = tuple(REGION_RATES[i * len(REGION_RATES) // length] for i in range(length))
    mutation_probabilities = [-expm1(-divergence * rate) for rate in site_rates]
    sequences = []
    for _ in range(size):
        residues = []
        for ancestral_residue, probability in zip(ancestor, mutation_probabilities):
            # Draw for each site at every divergence setting so the same seed
            # preserves ancestry and candidate changes across the sweep.
            draw = rng.random()
            alternative = rng.randrange(alphabet_size - 1)
            if alternative >= ancestral_residue:
                alternative += 1
            index = alternative if draw < probability else ancestral_residue
            residues.append(AMINO_ACIDS[index])
        sequences.append("".join(residues))

    ids = tuple(f"seq{i + 1}" for i in range(size))
    return SyntheticAlignment(
        alignment=Alignment(ids, tuple(sequences)),
        site_rates=site_rates,
        variable_sites=tuple(rate == REGION_RATES[-1] for rate in site_rates),
    )


def generate_constraint_alignment(
    *, size: int, length: int, divergence: float, seed: int
) -> SyntheticAlignment:
    """Pair within-class and cross-class substitutions at identical frequencies."""
    if size < 2 or length < 2 or length % 2:
        raise ValueError("the constraint suite needs at least two sequences and an even length")
    _check_divergence(divergence)
    rng = Random(seed)
    probability = -expm1(-divergence * MATCHED_SITE_RATE)
    rows = [[] for _ in range(size)]
    labels = []
    for _ in range(length // 2):
        group = rng.choice(RESIDUE_CLASSES)
        ancestral, preserved = rng.sample(group, 2)
        broken = rng.choice("".join(other for other in RESIDUE_CLASSES if other != group))
        crossed_first = rng.choice((False, True))
        replacements = (broken, preserved) if crossed_first else (preserved, broken)
        labels.extend((crossed_first, not crossed_first))
        for row in rows:
            pair = replacements if rng.random() < probability else (ancestral, ancestral)
            row.extend(pair)

    ids = tuple(f"seq{i + 1}" for i in range(size))
    return SyntheticAlignment(
        alignment=Alignment(ids, tuple("".join(row) for row in rows)),
        site_rates=(MATCHED_SITE_RATE,) * length,
        variable_sites=tuple(labels),
    )


def _check_divergence(divergence: float) -> None:
    if not isfinite(divergence) or divergence < 0:
        raise ValueError("divergence must be finite and nonnegative")

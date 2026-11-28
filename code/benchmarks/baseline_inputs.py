"""Provide INDELible rate inputs and matched-pair controls for the baseline benchmark."""

from dataclasses import dataclass
from math import expm1, isfinite
from pathlib import Path
from random import Random

from msa_conservation.msa import Alignment


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
    *, size: int, length: int, divergence: float, seed: int, work_dir: Path | None = None
) -> SyntheticAlignment:
    """Use INDELible for the original three-block rate comparison."""
    if size < 2 or length < 3:
        raise ValueError("synthetic alignments need at least two sequences and three sites")
    from data_generation.rvne_simulation import simulate, write_simulation

    rates = tuple(REGION_RATES[i * len(REGION_RATES) // length] for i in range(length))
    simulation = simulate(size=size, length=length, divergence=divergence, seed=seed,
                          ancestral_rates=rates, work_dir=work_dir)
    if work_dir is not None:
        write_simulation(simulation, work_dir / "aligned.fasta", work_dir / "truth.tsv", work_dir)
    return SyntheticAlignment(simulation.alignment, tuple(simulation.rates),
                              tuple(simulation.rates == REGION_RATES[-1]))


def generate_constraint_alignment(
    *, size: int, length: int, divergence: float, seed: int
) -> SyntheticAlignment:
    """Construct matched controls with equal mutation counts for both residue classes."""
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

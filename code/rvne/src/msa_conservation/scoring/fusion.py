"""Baseline entropy, substitution, and entropy-adjusted variability scores."""

from dataclasses import dataclass
from math import isfinite
from typing import Iterator, Sequence

from ..msa import Alignment, column_profiles
from .entropy import shannon_entropy
from .substitution import Blosum62, load_blosum62, substitution_scores


@dataclass(frozen=True)
class SiteScore:
    position: int
    nongap_coverage: float
    canonical_coverage: float
    entropy_variability: float | None
    substitution_variability: float | None
    conditional_severity: float | None
    integrated_variability: float | None


def score_alignment(
    alignment: Alignment,
    *,
    weights: Sequence[float] | None = None,
    gamma: float = 1.0,
    gap_penalty: float = 0.0,
    matrix: Blosum62 | None = None,
) -> Iterator[SiteScore]:
    """Yield site scores in alignment order. Higher means more variable."""
    if not isfinite(gamma) or gamma < 0:
        raise ValueError("gamma must be finite and nonnegative")
    if not isfinite(gap_penalty) or gap_penalty < 0:
        raise ValueError("gap penalty must be finite and nonnegative")
    if matrix is None:
        matrix = load_blosum62()
    for position, column in enumerate(column_profiles(alignment, weights), start=1):
        if column.canonical_coverage == 0:
            yield SiteScore(position, column.nongap_coverage, 0.0, None, None, None, None)
            continue

        entropy = shannon_entropy(column.frequencies)
        expected_cost, conditional = substitution_scores(column.frequencies, matrix)
        penalty = gap_penalty * (1 - column.nongap_coverage)
        yield SiteScore(
            position=position,
            nongap_coverage=column.nongap_coverage,
            canonical_coverage=column.canonical_coverage,
            entropy_variability=entropy + penalty,
            substitution_variability=expected_cost + penalty,
            conditional_severity=conditional,
            integrated_variability=entropy * (1 + gamma * conditional) + penalty,
        )

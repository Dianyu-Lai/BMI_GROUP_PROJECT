"""Published column-wise variability scores."""

from dataclasses import dataclass
from math import isfinite
from typing import Iterator, Sequence

from ..msa import Alignment, column_profiles
from .entropy import adjusted_entropy, jensen_shannon_divergence, relative_entropy, shannon_entropy
from .substitution import blosum62_sum_of_pairs, load_blosum62


@dataclass(frozen=True)
class SiteScore:
    position: int
    nongap_coverage: float
    canonical_coverage: float
    entropy_variability: float | None


@dataclass(frozen=True)
class ReferenceSiteScore(SiteScore):
    kl_variability: float | None
    jsd_variability: float | None
    blosum62_sop_variability: float | None
    adjusted_entropy_variability: float | None


def score_alignment(
    alignment: Alignment,
    *,
    weights: Sequence[float] | None = None,
    gap_penalty: float = 0.0,
    reference_baselines: bool = False,
) -> Iterator[SiteScore]:
    """Yield site scores in alignment order. Higher means more variable."""
    if not isfinite(gap_penalty) or gap_penalty < 0:
        raise ValueError("gap penalty must be finite and nonnegative")
    background = load_blosum62().background if reference_baselines else None
    profiles = column_profiles(alignment, weights, include_pair_mass=reference_baselines)
    for position, column in enumerate(profiles, start=1):
        if column.canonical_coverage == 0:
            missing = (position, column.nongap_coverage, 0.0, None)
            if reference_baselines:
                yield ReferenceSiteScore(*missing, None, None, None, None)
            else:
                yield SiteScore(*missing)
            continue

        frequencies = column.frequencies
        entropy = shannon_entropy(frequencies)
        penalty = gap_penalty * (1 - column.nongap_coverage)
        scores = (position, column.nongap_coverage, column.canonical_coverage, entropy + penalty)
        if background is None or column.same_residue_pair_mass is None:
            yield SiteScore(*scores)
            continue
        sop = blosum62_sum_of_pairs(frequencies, column.same_residue_pair_mass)
        yield ReferenceSiteScore(
            *scores,
            kl_variability=-relative_entropy(frequencies, background) + penalty,
            jsd_variability=-jensen_shannon_divergence(frequencies, background) + penalty,
            blosum62_sop_variability=-sop + penalty if sop is not None else None,
            adjusted_entropy_variability=adjusted_entropy(frequencies) + penalty,
        )

"""Coverage-weighted segmentation and adaptive smoothing for RVNE profiles."""

from dataclasses import dataclass
from math import isfinite

import numpy as np


@dataclass(frozen=True)
class ContextParameters:
    """Penalized segmentation and smoothing settings."""

    beta: float = 6.0
    smoothing: float = 0.25
    radii: tuple[int, ...] = (0, 1, 2, 4, 8)
    min_support: float = 1.5
    min_segment: int = 3
    max_segment: int | None = None

    def __post_init__(self) -> None:
        if not all(isfinite(x) and x >= 0 for x in (self.beta, self.smoothing, self.min_support)):
            raise ValueError("beta, smoothing and minimum support must be finite and nonnegative")
        if not self.radii or tuple(sorted(set(self.radii))) != self.radii or self.radii[0] != 0:
            raise ValueError("radii must be distinct and increasing, beginning at zero")
        if self.min_segment < 1 or (self.max_segment is not None and self.max_segment < self.min_segment):
            raise ValueError("invalid segment length constraints")


@dataclass(frozen=True)
class Segment:
    """A zero-based, half-open segment of the conservation profile."""

    start: int
    stop: int
    mean: float
    variance: float
    valid_columns: int


@dataclass(frozen=True)
class ContextScores:
    """Raw local estimates, coverage-penalized scores and per-site window metadata."""

    estimates: np.ndarray
    scores: np.ndarray
    window_size: np.ndarray
    variance: np.ndarray
    support: np.ndarray
    segment_id: np.ndarray
    segments: tuple[Segment, ...]


class WeightedIntervals:
    """Weighted interval statistics from centered prefix sums."""

    def __init__(self, scores: np.ndarray, coverage: np.ndarray) -> None:
        self.valid = np.isfinite(scores) & (coverage > 0)
        self.weights = np.where(self.valid, coverage, 0)
        values = np.where(self.valid, scores, 0)
        total = self.weights.sum()
        self.center = float(self.weights @ values / total) if total else 0.0
        centered = values - self.center
        self.weight = np.r_[0., np.cumsum(self.weights)]
        self.first = np.r_[0., np.cumsum(self.weights * centered)]
        self.second = np.r_[0., np.cumsum(self.weights * centered ** 2)]
        self.count = np.r_[0, np.cumsum(self.valid)]

    def moments(self, start, stop):
        """Return support, centered mean and weighted SSE on half-open intervals."""
        support = self.weight[stop] - self.weight[start]
        first = self.first[stop] - self.first[start]
        second = self.second[stop] - self.second[start]
        mean = first / np.where(support > 0, support, 1)
        return support, mean, np.maximum(second - first * mean, 0)

    @property
    def variance(self) -> float:
        return float(self.second[-1] / self.weight[-1]) if self.weight[-1] else 0.0

    def segment(self, start: int, stop: int) -> Segment:
        support, mean, error = self.moments(start, stop)
        return Segment(start, stop, float(self.center + mean), float(error / support),
                       int(self.count[stop] - self.count[start]))


def segment_profile(
    intervals: WeightedIntervals, params: ContextParameters
) -> tuple[Segment, ...]:
    """Solve penalized weighted least squares by optimal partitioning."""
    length = len(intervals.valid)
    if not intervals.count[-1]:
        return ()
    minimum = min(length, params.min_segment)
    maximum = params.max_segment or length
    penalty = params.beta * intervals.variance
    costs = np.full(length + 1, np.inf)
    previous = np.full(length + 1, -1, dtype=int)
    counts = np.zeros(length + 1, dtype=int)
    costs[0] = 0
    for stop in range(minimum, length + 1):
        starts = np.arange(max(0, stop - maximum), stop - minimum + 1)
        support, _, error = intervals.moments(starts, stop)
        objectives = np.where(support > 0, costs[starts] + error + penalty, np.inf)
        best = objectives.min()
        if not np.isfinite(best):
            continue
        tied = starts[objectives == best]
        start = tied[np.argmin(counts[tied])]
        costs[stop] = best
        previous[stop] = start
        counts[stop] = counts[start] + 1
    if previous[length] < 0:
        raise ValueError("no segmentation with positive support satisfies the length constraints")
    segments = []
    stop = length
    while stop:
        start = int(previous[stop])
        segments.append(intervals.segment(start, stop))
        stop = start
    return tuple(reversed(segments))


METHODS = ("uniform", "gaussian", "adaptive", "segmentation", "scaw")


def contextual_scores(
    raw: np.ndarray, coverage: np.ndarray, method: str,
    params: ContextParameters = ContextParameters(),
) -> ContextScores:
    """Score a profile using one context method, preserving missing centers."""
    if method not in METHODS:
        raise ValueError(f"unknown context method: {method}")
    intervals = WeightedIntervals(raw, coverage)
    length = len(raw)
    estimates = np.full(length, np.nan)
    variance = estimates.copy()
    support = np.zeros(length)
    width = np.zeros(length, dtype=int)
    segment_id = width.copy()
    segments = segment_profile(intervals, params) if method in ("segmentation", "scaw") else ()
    if not intervals.count[-1]:
        return ContextScores(estimates, estimates.copy(), width, variance, support, segment_id, segments)
    regions = segments or (intervals.segment(0, length),)
    radii = np.asarray(params.radii)
    radius = radii[-1]
    scale = intervals.variance
    for label, segment in enumerate(regions, 1):
        if segments:
            segment_id[segment.start:segment.stop] = label
        positions = np.flatnonzero(intervals.valid[segment.start:segment.stop]) + segment.start
        if method == "segmentation":
            estimates[positions] = segment.mean
            variance[positions] = segment.variance
            support[positions] = intervals.weight[segment.stop] - intervals.weight[segment.start]
            width[positions] = segment.stop - segment.start
        elif method in ("adaptive", "scaw"):
            starts = np.maximum(segment.start, positions[:, None] - radii)
            stops = np.minimum(segment.stop, positions[:, None] + radii + 1)
            weights, means, errors = intervals.moments(starts, stops)
            variances = errors / weights
            objective = params.smoothing / weights
            if scale > 0:
                objective += variances / scale
            objective[:, 1:] = np.where(weights[:, 1:] >= params.min_support,
                                        objective[:, 1:], np.inf)
            # Search widest first, retaining the largest candidate on exact ties.
            selected = len(radii) - 1 - np.argmin(objective[:, ::-1], axis=1)
            rows = np.arange(len(positions))
            pick = rows, selected
            estimates[positions] = intervals.center + means[pick]
            variance[positions] = variances[pick]
            support[positions] = weights[pick]
            width[positions] = stops[pick] - starts[pick]
        elif method == "uniform":
            starts = np.maximum(positions - radius, 0)
            stops = np.minimum(positions + radius + 1, length)
            weights, means, errors = intervals.moments(starts, stops)
            estimates[positions] = intervals.center + means
            variance[positions] = errors / weights
            support[positions] = weights
            width[positions] = stops - starts
        else:
            offsets = np.arange(-radius, radius + 1)
            kernel = np.exp(-0.5 * (offsets / max(0.5, radius / 2)) ** 2)
            for i in positions:
                start, stop = max(0, i - radius), min(i + radius + 1, length)
                weights = intervals.weights[start:stop] * kernel[start - i + radius:stop - i + radius]
                values = np.nan_to_num(raw[start:stop], nan=0.0)
                support[i] = weights.sum()
                estimates[i] = weights @ values / support[i]
                variance[i] = weights @ ((values - estimates[i]) ** 2) / support[i]
                width[i] = stop - start
    return ContextScores(estimates, estimates * coverage, width, variance, support, segment_id, segments)

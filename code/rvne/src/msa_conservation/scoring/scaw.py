"""SCAW v2: penalized weighted SSE and a fixed, boundary-clipped context window."""

from dataclasses import dataclass
from math import isfinite

import numpy as np


@dataclass(frozen=True)
class ContextParameters:
    """One multiplier of the data-derived BIC-style penalty (not old beta)."""

    penalty_scale: float = 1.0

    def __post_init__(self) -> None:
        if (type(self.penalty_scale) not in (int, float) or
                not isfinite(self.penalty_scale) or self.penalty_scale < 0):
            raise ValueError("penalty_scale must be finite and nonnegative")


METHOD_VERSION = "scaw-v2"
WINDOW_RADIUS = 3


def context_from_settings(values: dict) -> ContextParameters:
    """Reject historical settings rather than silently changing their meaning."""
    if not isinstance(values, dict) or set(values) - {"penalty_scale"}:
        raise ValueError("SCAW v2 accepts only penalty_scale; v1 settings are unsupported")
    return ContextParameters(**values)


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
        scores, coverage = np.asarray(scores, dtype=float), np.asarray(coverage, dtype=float)
        if scores.ndim != 1 or scores.shape != coverage.shape:
            raise ValueError("scores and coverage must be equal-length one-dimensional arrays")
        if not np.all(np.isfinite(coverage) & (coverage >= 0) & (coverage <= 1)):
            raise ValueError("coverage must be finite and between zero and one")
        if np.isinf(scores).any():
            raise ValueError("scores may contain NaN, but not infinity")
        self.scores = scores
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


def noise_variance(intervals: WeightedIntervals) -> float:
    """First-difference variance, using only adjacent scorable original columns.

    Under independent errors Var(y_i)=sigma^2/w_i, E[(y_i-y_j)^2 /
    (1/w_i+1/w_j)] = sigma^2 within a constant regime. For unit weights this
    reduces to Rice's first-difference estimator. Changes bias it upward;
    positive serial correlation biases it downward. See the methods document.
    """
    adjacent = intervals.valid[:-1] & intervals.valid[1:]
    if not adjacent.any():
        return intervals.variance
    difference = np.diff(intervals.scores)[adjacent]
    weights = intervals.weights
    denominator = 1 / weights[:-1][adjacent] + 1 / weights[1:][adjacent]
    return float(np.mean(difference ** 2 / denominator))


def segmentation_penalty(intervals: WeightedIntervals, params: ContextParameters) -> float:
    """Return c * sigma_hat^2 * 2 log(m), in weighted SSE units."""
    return float(params.penalty_scale * noise_variance(intervals) *
                 (2 * np.log(max(2, int(intervals.count[-1])))))


def segment_profile(
    intervals: WeightedIntervals, params: ContextParameters
) -> tuple[Segment, ...]:
    """Exact PELT for sum(weighted SSE) + penalty * number of changes.

    Optimize on scorable columns; map changes to the next scorable column.
    Missing runs carry no evidence for boundaries and stay with the left segment
    (leading missing columns belong to the first segment). No length constraints.
    """
    length = len(intervals.valid)
    if not intervals.count[-1]:
        return ()
    indices = np.flatnonzero(intervals.valid)
    compact = WeightedIntervals(intervals.scores[indices], intervals.weights[indices])
    size = len(indices)
    penalty = segmentation_penalty(intervals, params)
    costs = np.full(size + 1, np.inf)
    previous = np.full(size + 1, -1, dtype=int)
    counts = np.zeros(size + 1, dtype=int)
    costs[0] = -penalty
    candidates = np.array([0], dtype=int)
    # Prefix-sum roundoff tolerance, not an inferential tuning parameter.
    tolerance = 64 * np.finfo(float).eps * max(float(compact.second[-1]), penalty,
                                               np.finfo(float).tiny)
    for stop in range(1, size + 1):
        starts = candidates
        _, _, error = compact.moments(starts, stop)
        unpenalized = costs[starts] + error
        objectives = unpenalized + penalty
        best = objectives.min()
        if not np.isfinite(best):
            continue
        tied = starts[objectives <= best + tolerance]
        start = tied[np.argmin(counts[tied])]
        costs[stop] = objectives[np.flatnonzero(starts == start)[0]]
        previous[stop] = start
        counts[stop] = counts[start] + 1
        # Weighted SSE is superadditive, so PELT's pruning constant is zero.
        candidates = np.r_[starts[unpenalized <= costs[stop] + tolerance], stop]
    segments = []
    stop = size
    while stop:
        start = int(previous[stop])
        left = int(indices[start]) if start else 0
        right = int(indices[stop]) if stop < size else length
        segments.append(intervals.segment(left, right))
        stop = start
    return tuple(reversed(segments))


METHODS = ("uniform", "gaussian", "segmentation", "scaw")


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
    radius = WINDOW_RADIUS
    for label, segment in enumerate(regions, 1):
        if segments:
            segment_id[segment.start:segment.stop] = label
        positions = np.flatnonzero(intervals.valid[segment.start:segment.stop]) + segment.start
        if method == "segmentation":
            estimates[positions] = segment.mean
            variance[positions] = segment.variance
            support[positions] = intervals.weight[segment.stop] - intervals.weight[segment.start]
            width[positions] = segment.stop - segment.start
        elif method in ("uniform", "scaw"):
            starts = np.maximum(positions - radius, segment.start)
            stops = np.minimum(positions + radius + 1, segment.stop)
            weights, means, errors = intervals.moments(starts, stops)
            estimates[positions] = intervals.center + means
            if method == "scaw":
                # Fixed 1/2 center + 1/2 coverage-weighted neighbors, following
                # Capra & Singh's window heuristic with our boundary/mask rules.
                neighbor_weight = weights - intervals.weights[positions]
                neighbor_first = (means * weights - intervals.weights[positions] *
                                  (raw[positions] - intervals.center))
                neighbor_mean = intervals.center + np.divide(
                    neighbor_first, neighbor_weight, out=np.zeros_like(neighbor_weight),
                    where=neighbor_weight > 0)
                estimates[positions] = np.where(
                    neighbor_weight > 0, .5 * (raw[positions] + neighbor_mean), raw[positions])
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

"""Evaluation shared by the RVNE parameter search and held-out benchmark."""

from itertools import groupby

import numpy as np

from msa_conservation.msa import AMINO_ACIDS, Alignment, column_profiles
from msa_conservation.scoring.rvne import LOG_SIZE, SiteProfiles, background
from msa_conservation.scoring.substitution import load_blosum62
from .metrics import roc_auc


def _ranks(values: np.ndarray) -> np.ndarray:
    """Vectorized average ranks along the final axis, including tied values."""
    order = np.argsort(values, axis=-1, kind="stable")
    sorted_values = np.take_along_axis(values, order, axis=-1)
    indices = np.arange(values.shape[-1])
    different = sorted_values[..., 1:] != sorted_values[..., :-1]
    starts = np.maximum.accumulate(np.where(
        np.concatenate((np.ones_like(sorted_values[..., :1], dtype=bool), different), axis=-1),
        indices, 0), axis=-1)
    ends = np.minimum.accumulate(np.where(
        np.concatenate((different, np.ones_like(sorted_values[..., :1], dtype=bool)), axis=-1),
        indices, len(indices) - 1)[..., ::-1], axis=-1)[..., ::-1]
    ranks = np.empty_like(values, dtype=float)
    np.put_along_axis(ranks, order, (starts + ends) / 2, axis=-1)
    return ranks


def correlations(values: np.ndarray, truth: np.ndarray, *, rank: bool = False) -> np.ndarray:
    """Correlate each row with one shared target, allowing undefined rows."""
    if truth.size < 2:
        return np.full(values.shape[:-1], np.nan)
    if rank:
        values, truth = _ranks(np.round(values, 12)), _ranks(truth)
    x = values - values.mean(axis=-1, keepdims=True)
    y = truth - truth.mean()
    denom = np.linalg.norm(x, axis=-1) * np.linalg.norm(y)
    return np.divide(x @ y, denom, out=np.full(np.shape(denom), np.nan), where=denom > 1e-15)


def average_precision(scores: np.ndarray, labels: np.ndarray) -> float:
    """Compute step-integrated precision-recall area, grouping equal scores."""
    if not np.any(labels) or np.all(labels):
        return float("nan")
    tp = seen = 0
    area = 0.0
    pairs = sorted(zip(scores, labels), key=lambda pair: pair[0], reverse=True)
    for _, group in groupby(pairs, key=lambda pair: pair[0]):
        group_labels = [bool(label) for _, label in group]
        cases = sum(group_labels)
        tp += cases
        seen += len(group_labels)
        area += cases * tp / seen
    return area / int(labels.sum())


def site_metrics(scores: np.ndarray, truth: np.ndarray, labels: np.ndarray,
                 threshold: float) -> dict[str, float]:
    """Measure site performance on the scores the method actually emitted."""
    valid = np.isfinite(scores)
    scores = np.round(scores[valid], 12)
    truth, labels = truth[valid], labels[valid]
    positive, negative = int(labels.sum()), int((~labels).sum())
    enough = positive and negative
    return {
        "scored_fraction": float(valid.mean()),
        "spearman": float(correlations(scores, truth, rank=True)),
        "pearson": float(correlations(scores, truth)),
        "auroc": roc_auc(scores.tolist(), labels.tolist()) if enough else float("nan"),
        "auprc": average_precision(scores, labels),
        "sensitivity": float(np.count_nonzero(labels & (scores >= threshold)) / positive)
        if positive else float("nan"),
        "specificity": float(np.count_nonzero(~labels & (scores < threshold)) / negative)
        if negative else float("nan"),
    }


def boundary_metrics(predicted: tuple[int, ...], truth: tuple[int, ...],
                     tolerance: int = 2) -> tuple[float, float, float]:
    """Match sorted boundaries one-to-one within a prespecified tolerance."""
    i = j = matches = 0
    while i < len(predicted) and j < len(truth):
        delta = predicted[i] - truth[j]
        if abs(delta) <= tolerance:
            matches += 1
            i += 1
            j += 1
        elif delta < 0:
            i += 1
        else:
            j += 1
    if not predicted and not truth:
        return 1., 1., 1.
    return (matches / len(predicted) if predicted else 0.,
            matches / len(truth) if truth else 0.,
            2 * matches / (len(predicted) + len(truth)))


def region_iou(prediction: np.ndarray, labels: np.ndarray, valid: np.ndarray) -> float:
    """Intersection over union of contiguous positive masks on evaluated coordinates."""
    union = np.count_nonzero((prediction | labels) & valid)
    return float(np.count_nonzero(prediction & labels & valid) / union) if union else float("nan")


# ClustalW conservation symbols, applied to the supplied alignment without realigning it.
_CLUSTAL_STRONG = tuple(map(frozenset, "STA NEQK NHQK NDEQ QHRK MILV MILF HY FYW".split()))
_CLUSTAL_WEAK = tuple(map(frozenset, "CSA ATV SAG STNK STPA SGND SNDEQK NDEQHK NEQHRK FVLIM HFY".split()))


def _clustal_consensus(column: tuple[str, ...]) -> int:
    residues = set(column)
    if not residues.issubset(AMINO_ACIDS):
        return 0
    if len(residues) == 1:
        return 3
    if any(residues <= group for group in _CLUSTAL_STRONG):
        return 2
    return int(any(residues <= group for group in _CLUSTAL_WEAK))


def _al2co_variance(alignment: Alignment) -> np.ndarray:
    """Score the unweighted AL2CO variance measure using pooled canonical counts."""
    profiles = list(column_profiles(alignment))
    probabilities = np.array([site.frequencies for site in profiles])
    canonical_coverage = np.array([site.canonical_coverage for site in profiles])
    total = canonical_coverage.sum()
    if total == 0:
        return np.full(alignment.length, np.nan)
    pooled = (canonical_coverage @ probabilities) / total
    return np.square(probabilities - pooled).sum(axis=1)


def comparator_profiles(sites: SiteProfiles, alignment: Alignment) -> dict[str, np.ndarray]:
    """Score column comparators with the shared RVNE mask and coverage multiplier."""
    p = sites.probabilities
    q = background(1)
    middle = (p + q) / 2
    positive = p > 0
    log_p = np.log(p, where=positive, out=np.zeros_like(p))
    shannon = 1 + (p * log_p).sum(axis=1) / LOG_SIZE
    kl = (p * (log_p - np.log(q))).sum(axis=1) / np.log(2)
    jsd = ((p * (log_p - np.log(middle))).sum(axis=1)
           + (q * np.log(q / middle)).sum(axis=1)) / (2 * np.log(2))
    blosum = load_blosum62()
    severity = np.einsum("li,ij,lj->l", p, np.asarray(blosum.severity), p)
    similarity = 1 - severity / np.max(blosum.severity)
    clustal = np.fromiter(map(_clustal_consensus, zip(*alignment.sequences)),
                          dtype=float, count=alignment.length)
    return {method: np.where(sites.valid, values * sites.coverage, np.nan)
            for method, values in (("Shannon", shannon), ("JSD", jsd), ("BLOSUM62", similarity),
                                   ("KL", kl), ("ClustalW consensus", clustal),
                                   ("AL2CO variance", _al2co_variance(alignment)))}

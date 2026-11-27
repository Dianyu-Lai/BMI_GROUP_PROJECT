"""Tie-aware ranking and binary classification metrics for synthetic site labels."""

from itertools import groupby
from math import fsum, inf, sqrt
from typing import Sequence


def _average_ranks(values: Sequence[float]) -> list[float]:
    indices = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    for _, group in groupby(indices, key=values.__getitem__):
        tied = list(group)
        end = start + len(tied)
        average_rank = (start + 1 + end) / 2
        for i in tied:
            ranks[i] = average_rank
        start = end
    return ranks


def spearman(scores: Sequence[float], site_rates: Sequence[float]) -> float | None:
    if len(scores) != len(site_rates) or not scores:
        raise ValueError("scores and rates must have the same nonzero length")
    center = (len(scores) + 1) / 2
    x = [rank - center for rank in _average_ranks(scores)]
    y = [rank - center for rank in _average_ranks(site_rates)]
    scale = sqrt(fsum(v * v for v in x) * fsum(v * v for v in y))
    return fsum(a * b for a, b in zip(x, y)) / scale if scale else None


def _class_counts(scores: Sequence[float], labels: Sequence[bool]) -> tuple[int, int]:
    if len(scores) != len(labels):
        raise ValueError("scores and binary labels must have equal length")
    positives = sum(labels)
    negatives = len(labels) - positives
    if not positives or not negatives:
        raise ValueError("binary labels must contain both classes")
    return positives, negatives


def roc_auc(scores: Sequence[float], labels: Sequence[bool]) -> float:
    positives, negatives = _class_counts(scores, labels)
    positive_rank_sum = fsum(rank for rank, label in zip(_average_ranks(scores), labels) if label)
    return (positive_rank_sum - positives * (positives + 1) / 2) / (positives * negatives)


def select_threshold(scores: Sequence[float], labels: Sequence[bool]) -> float:
    """Maximize balanced accuracy, preferring the higher cutoff on ties."""
    positives, negatives = _class_counts(scores, labels)
    ranked = sorted(zip(scores, labels), key=lambda sample: sample[0], reverse=True)
    best_threshold = inf
    best_objective = 0
    tp = fp = 0
    for score, group in groupby(ranked, key=lambda sample: sample[0]):
        for _, label in group:
            tp += label
            fp += not label
        # Maximizing TP / P - FP / N also maximizes balanced accuracy.
        objective = tp * negatives - fp * positives
        if objective > best_objective:
            best_objective = objective
            best_threshold = score
    return best_threshold


def sensitivity_specificity(
    scores: Sequence[float], labels: Sequence[bool], threshold: float
) -> tuple[float, float]:
    positives, negatives = _class_counts(scores, labels)
    tp = sum(score >= threshold and label for score, label in zip(scores, labels))
    tn = sum(score < threshold and not label for score, label in zip(scores, labels))
    return tp / positives, tn / negatives


def matched_pair_accuracy(scores: Sequence[float], labels: Sequence[bool]) -> float:
    """Score cross-class columns above their matched controls, with half credit for ties."""
    if len(scores) != len(labels) or not scores or len(scores) % 2:
        raise ValueError("matched scores and labels must have the same positive even length")
    wins = 0.0
    for i in range(0, len(scores), 2):
        if labels[i] == labels[i + 1]:
            raise ValueError("each matched pair must contain both classes")
        positive, negative = (scores[i], scores[i + 1]) if labels[i] else (scores[i + 1], scores[i])
        if positive > negative:
            wins += 1.0
        elif positive == negative:
            wins += 0.5
    return wins / (len(scores) // 2)

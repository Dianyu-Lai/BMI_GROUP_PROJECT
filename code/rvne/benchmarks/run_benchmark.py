"""Score separate calibration and test alignments in two controlled synthetic suites."""

import argparse
import csv
import sys
from time import perf_counter

from msa_conservation.scoring.fusion import score_alignment
from msa_conservation.scoring.physicochemical import calibrated_blosum62
from msa_conservation.scoring.substitution import Blosum62, load_blosum62

from .metrics import matched_pair_accuracy, roc_auc, select_threshold, sensitivity_specificity, spearman
from .synthetic import SyntheticAlignment, generate_alignment, generate_constraint_alignment


METHODS = ("entropy_variability", "substitution_variability", "integrated_variability")


def _score_columns(
    dataset: SyntheticAlignment,
    gamma: float,
    matrix: Blosum62 | None = None,
    methods: tuple[str, ...] = METHODS,
) -> tuple[list[float], ...]:
    columns = tuple([] for _ in methods)
    for site in score_alignment(dataset.alignment, gamma=gamma, matrix=matrix):
        for column, method in zip(columns, methods):
            value = getattr(site, method)
            if value is None:
                raise ValueError(f"synthetic site {site.position} has no score")
            column.append(value)
    return columns


def _benchmark_scores(
    dataset: SyntheticAlignment, gamma: float, matrix: Blosum62 | None
) -> tuple[list[float], ...]:
    scores = _score_columns(dataset, gamma)
    if matrix is not None:
        scores += _score_columns(dataset, gamma, matrix, METHODS[1:])
    return scores


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=64, help="number of sequences per alignment")
    parser.add_argument("--length", type=int, default=300, help="number of aligned sites")
    parser.add_argument("--divergence", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=1, help="calibration seed (test uses seed + 1)")
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--suite", choices=("rate", "constraint"), default="rate")
    parser.add_argument(
        "--aaindex-eta", type=float,
        help="include AAindex-calibrated scores (default: 0.5 for constraint, off for rate)",
    )
    args = parser.parse_args()

    parameters = dict(size=args.size, length=args.length, divergence=args.divergence)
    constraint = args.suite == "constraint"
    generator = generate_constraint_alignment if constraint else generate_alignment
    eta = args.aaindex_eta
    if eta is None and constraint:
        eta = 0.5
    try:
        calibration = generator(**parameters, seed=args.seed)
        test = generator(**parameters, seed=args.seed + 1)
        load_blosum62()
        matrix = calibrated_blosum62(eta) if eta is not None else None
        calibration_scores = _benchmark_scores(calibration, args.gamma, matrix)
        start = perf_counter()
        test_scores = _benchmark_scores(test, args.gamma, matrix)
        elapsed = perf_counter() - start
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    header = ("method", "spearman", "auroc", "sensitivity", "specificity", "threshold")
    if constraint:
        header += ("matched_pair_accuracy",)
    writer.writerow(header)
    methods = METHODS
    if matrix is not None:
        methods += ("aaindex_substitution_variability", "aaindex_integrated_variability")
    truth = test.variable_sites if constraint else test.site_rates
    for method, fitted, observed in zip(methods, calibration_scores, test_scores):
        threshold = select_threshold(fitted, calibration.variable_sites)
        sens, spec = sensitivity_specificity(observed, test.variable_sites, threshold)
        row = (
            method, spearman(observed, truth), roc_auc(observed, test.variable_sites),
            sens, spec, threshold,
        )
        if constraint:
            row += (matched_pair_accuracy(observed, test.variable_sites),)
        writer.writerow(row)
    print(f"Test scoring time for {len(methods)} methods: {elapsed:.4f} s", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

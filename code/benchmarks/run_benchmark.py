"""Score separate calibration and test alignments in two controlled synthetic suites."""

import argparse
import csv
import sys
from time import perf_counter
from pathlib import Path

from .baseline_inputs import (
    SyntheticAlignment, generate_alignment, generate_constraint_alignment,
)
from msa_conservation.scoring.baseline import score_alignment
from msa_conservation.scoring.substitution import load_blosum62

from .metrics import matched_pair_accuracy, roc_auc, select_threshold, sensitivity_specificity, spearman


METHODS = ("entropy_variability",)
REFERENCE_METHODS = ("kl_variability", "jsd_variability", "blosum62_sop_variability", "adjusted_entropy_variability")


def _score_columns(
    dataset: SyntheticAlignment,
    methods: tuple[str, ...],
    *,
    include_references: bool,
) -> tuple[list[float], ...]:
    columns = tuple([] for _ in methods)
    for site in score_alignment(dataset.alignment, reference_baselines=include_references):
        for column, method in zip(columns, methods):
            value = getattr(site, method)
            if value is None:
                raise ValueError(f"synthetic site {site.position} has no score")
            column.append(value)
    return columns


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=64, help="number of sequences per alignment")
    parser.add_argument("--length", type=int, default=300, help="number of ancestral sites")
    parser.add_argument("--divergence", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=1, help="calibration seed (test uses seed + 1)")
    parser.add_argument("--suite", choices=("rate", "constraint"), default="rate")
    parser.add_argument("--data-output", type=Path, help="new directory for evolutionary inputs and metadata")
    parser.add_argument(
        "--reference-baselines", action="store_true", help="compare KL, JSD, BLOSUM62 SOP and adjusted entropy"
    )
    args = parser.parse_args()

    parameters = dict(size=args.size, length=args.length, divergence=args.divergence)
    constraint = args.suite == "constraint"
    if args.data_output and (constraint or args.data_output.exists()):
        parser.error("--data-output must be a new directory for the INDELible rate suite")
    generator = generate_constraint_alignment if constraint else generate_alignment
    methods = METHODS + REFERENCE_METHODS if args.reference_baselines else METHODS
    try:
        calibration_output = {"work_dir": args.data_output / "calibration"} if args.data_output else {}
        test_output = {"work_dir": args.data_output / "test"} if args.data_output else {}
        calibration = generator(**parameters, seed=args.seed, **calibration_output)
        test = generator(**parameters, seed=args.seed + 1, **test_output)
        if args.reference_baselines:
            load_blosum62()
        calibration_scores = _score_columns(calibration, methods, include_references=args.reference_baselines)
        start = perf_counter()
        test_scores = _score_columns(test, methods, include_references=args.reference_baselines)
        elapsed = perf_counter() - start
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    header = ("method", "spearman", "auroc", "sensitivity", "specificity", "threshold")
    if constraint:
        header += ("matched_pair_accuracy",)
    writer.writerow(header)
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

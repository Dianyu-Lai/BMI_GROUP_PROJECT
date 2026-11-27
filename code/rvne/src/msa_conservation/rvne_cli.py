"""Score an aligned protein FASTA and export site and conservation-segment TSVs."""

import argparse
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys
from typing import Sequence

import numpy as np

from .msa import read_fasta
from .scoring.rvne import DEFAULT_PARAMETERS, RVNEParameters, score_reference, site_profiles
from .scoring.rvne_optimized import score_optimized
from .scoring.scaw import ContextParameters, METHODS, contextual_scores


def _number(value: float) -> str:
    return f"{value:.12g}" if np.isfinite(value) else ""


def _write_table(path: Path | None, header: tuple[str, ...], rows) -> None:
    if path is None:
        writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _site_rows(raw, sites, fitted, params):
    for i, score in enumerate(raw):
        valid = sites.valid[i]
        segment_id = int(fitted.segment_id[i]) if fitted is not None else 0
        segment = fitted.segments[segment_id - 1] if segment_id else None
        estimate = fitted.estimates[i] if fitted is not None else score
        yield (i + 1, _number(score), _number(sites.coverage[i]), _number(estimate),
               _number(estimate * sites.coverage[i]),
               int(fitted.window_size[i]) if fitted is not None and valid else (1 if valid else ""),
               _number(fitted.variance[i]) if fitted is not None else "",
               _number(fitted.support[i]) if fitted is not None and valid else "",
               segment_id or "", segment.start + 1 if segment else "", segment.stop if segment else "",
               params.r, params.s1, params.s2)


def _segment_rows(fitted, lower: float | None, upper: float | None):
    for i, segment in enumerate(fitted.segments, 1):
        label = "unclassified"
        if lower is not None and upper is not None:
            label = "variable" if segment.mean < lower else "intermediate"
            if segment.mean >= upper:
                label = "conserved"
        yield (i, segment.start + 1, segment.stop, _number(segment.mean),
               _number(segment.variance), segment.valid_columns, label)


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the score subcommand."""
    parser = argparse.ArgumentParser(prog="msa_conservation.cli score", description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="aligned protein FASTA")
    parser.add_argument("--output", type=Path, help="site TSV, or stdout when omitted")
    parser.add_argument("--segments-output", type=Path, help="optional segment TSV")
    parser.add_argument("--settings-output", type=Path, help="optional JSON with effective parameters")
    parser.add_argument("--method", choices=("rvne", "rvne-optimized", *METHODS[:-1], "scaw-rvne"),
                        default="rvne")
    parser.add_argument("--settings", type=Path, help="frozen parameter JSON from the benchmark")
    for name in ("r", "s1", "s2"):
        parser.add_argument(f"--{name}", type=float)
    parser.add_argument("--max-gap", type=float, default=.3)
    parser.add_argument("--beta", type=float)
    parser.add_argument("--smoothing", type=float, help="lambda in the window objective")
    parser.add_argument("--radii", type=int, nargs="+")
    parser.add_argument("--min-support", type=float)
    parser.add_argument("--min-segment", type=int)
    parser.add_argument("--max-segment", type=int, help="use zero for an unbounded partition")
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--segment-thresholds", type=float, nargs=2, metavar=("LOWER", "UPPER"))
    args = parser.parse_args(argv)
    lower, upper = args.segment_thresholds or (None, None)
    if lower is not None and (not np.isfinite([lower, upper]).all() or not lower < upper):
        parser.error("segment thresholds must be finite and strictly increasing")
    if (args.segments_output or args.segment_thresholds) and args.method not in ("segmentation", "scaw-rvne"):
        parser.error("segment export and thresholds require segmentation or scaw-rvne")
    paths = [path.resolve() for path in (args.input, args.settings, args.output,
                                       args.segments_output, args.settings_output) if path is not None]
    if len(set(paths)) != len(paths):
        parser.error("input, settings and output files must be distinct")
    try:
        settings = json.loads(args.settings.read_text(encoding="utf-8")) if args.settings else {}
        if not isinstance(settings, dict) or not isinstance(settings.get("context", {}), dict):
            raise ValueError("settings and context must be JSON objects")
        params = replace(RVNEParameters(**settings.get("rvne", asdict(DEFAULT_PARAMETERS))),
                         **{key: getattr(args, key) for key in ("r", "s1", "s2")
                            if getattr(args, key) is not None})
        context_values = settings.get("context", {}).copy()
        for key in ("beta", "smoothing", "radii", "min_support", "min_segment", "max_segment"):
            value = getattr(args, key)
            if value is not None:
                context_values[key] = None if key == "max_segment" and value == 0 else value
        if "radii" in context_values:
            context_values["radii"] = tuple(context_values["radii"])
        context = ContextParameters(**context_values)
        with args.input.open(encoding="utf-8-sig") as handle:
            sites = site_profiles(read_fasta(handle), args.max_gap)
        raw = (score_reference(sites, params) if args.method == "rvne" else
               score_optimized(sites, params, batch_size=args.batch_size))
        method = "scaw" if args.method == "scaw-rvne" else args.method
        fitted = contextual_scores(raw, sites.coverage, method, context) if method in METHODS else None
        _write_table(args.output,
                     ("position", "raw_rvne", "nongap_coverage", "local_estimate", "score", "window_size",
                      "local_variance", "window_support", "segment_id", "segment_start", "segment_end",
                      "r", "s1", "s2"),
                     _site_rows(raw, sites, fitted, params))
        if args.segments_output:
            _write_table(args.segments_output,
                         ("segment_id", "start", "end", "raw_mean", "variance", "valid_columns", "class"),
                         _segment_rows(fitted, lower, upper))
        if args.settings_output:
            effective = {"method": args.method, "rvne": asdict(params), "context": asdict(context),
                         "max_gap": args.max_gap, "batch_size": args.batch_size,
                         "segment_thresholds": args.segment_thresholds}
            args.settings_output.write_text(json.dumps(effective, indent=2) + "\n", encoding="utf-8")
    except (OSError, ValueError, TypeError, np.linalg.LinAlgError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Tune and evaluate RVNE variants on independent synthetic protein alignments."""

import argparse
import csv
import json
import platform
from collections import defaultdict
from dataclasses import asdict, dataclass, replace
from itertools import product
from pathlib import Path
from statistics import mean, median, stdev
from time import perf_counter
import tracemalloc

import numpy as np

from msa_conservation.scoring.rvne import (
    DEFAULT_PARAMETERS, RVNEParameters, SiteProfiles, paper_parameter_grid,
    score_reference, site_profiles,
)
from msa_conservation.scoring.rvne_optimized import score_grid, score_optimized
from msa_conservation.scoring.scaw import ContextParameters, contextual_scores
from .metrics import select_threshold
from .rvne_evaluation import boundary_metrics, comparator_profiles, correlations, region_iou, site_metrics
from .rvne_simulation import SCENARIOS, Simulation, simulate


CONTEXT = ContextParameters()
BASE_SPECS = tuple({"scenario": name} for name in SCENARIOS)
TRAINING = (101, 102, 103)
VALIDATION = (201, 202, 203)
HELD_OUT = (301, 302, 303, 304)
CONTEXT_METHODS = {"E2": "uniform", "E3": "gaussian", "E4": "adaptive",
                   "E5": "segmentation", "E6": "scaw"}
SWEEPS = (
    *({"size": n} for n in (20, 100, 500)),
    *({"length": length} for length in (100, 500, 1000, 5000)),
    *({"divergence": d} for d in (.25, 1.4)),
    *({"scenario": "motifs", "motif_length": n} for n in (3, 12)),
    *({"variable_fraction": p} for p in (.2, .6, .8)),
    *({"scenario": "indels", "gaps": p} for p in (.1, .3, .5)),
    *({"scenario": "redundant", "redundancy": p} for p in (.3, .85)),
    *({"mutation": mode} for mode in ("conservative", "radical")),
)


@dataclass(frozen=True)
class Family:
    seed: int
    scenario: str
    settings: dict
    truth: Simulation
    sites: SiteProfiles


def families(specs: tuple[dict, ...], seed_groups: tuple[int, ...]) -> list[Family]:
    """Generate independent alignments with unique RNG seeds within a split."""
    result = []
    for group, (index, spec) in product(seed_groups, enumerate(specs)):
        seed = group * 100 + index
        sim = simulate(seed=seed, **spec)
        name = spec.get("scenario", "blocks")
        result.append(Family(seed, name, spec, sim, site_profiles(sim.alignment)))
    return result


def grid_objectives(inputs: list[Family], grid: tuple[RVNEParameters, ...]) -> np.ndarray:
    """Return family-by-setting Spearman correlations, excluding uniform targets."""
    rows = []
    for family in inputs:
        select = family.sites.valid
        truth = family.truth.truth[select]
        if len(np.unique(truth)) > 1:
            scores = score_grid(family.sites, grid) * family.sites.coverage[select]
            rows.append(correlations(scores, truth, rank=True))
    if not rows:
        raise ValueError("parameter search needs families with variable ground truth and usable columns")
    return np.asarray(rows)


def _averages(values: np.ndarray) -> np.ndarray:
    available = np.isfinite(values).sum(axis=0)
    return np.divide(np.nansum(values, axis=0), available,
                     out=np.full(values.shape[1], -np.inf), where=available > 0)


def _sorted(values: np.ndarray, grid: tuple[RVNEParameters, ...]) -> list[int]:
    return sorted(range(len(grid)), key=lambda i: (
        -round(float(values[i]), 12), sum(abs(a-b) for a,b in zip(
            asdict(grid[i]).values(), asdict(DEFAULT_PARAMETERS).values())), grid[i]))


def select_parameters(training: list[Family], validation: list[Family]) -> tuple[dict, list[dict]]:
    """Select a global RVNE setting, then segment and window settings."""
    coarse = paper_parameter_grid()
    start = perf_counter()
    by_family = grid_objectives(training, coarse)
    grid_seconds = perf_counter() - start
    coarse_means = _averages(by_family)
    ranking = _sorted(coarse_means, coarse)
    best_training = coarse[ranking[0]]
    local_values = []
    for parameter, limit in zip(asdict(best_training).values(), (1., 1., .99)):
        local_values.append(sorted({round(min(limit, max(0., parameter + delta)), 4)
                                    for delta in (-.05, 0., .05)}))
    candidates = tuple(sorted({DEFAULT_PARAMETERS, *(coarse[i] for i in ranking[:12]),
                               *(RVNEParameters(*vals) for vals in product(*local_values))}))
    validation_objectives = grid_objectives(validation, candidates)
    validation_means = _averages(validation_objectives)
    order = _sorted(validation_means, candidates)
    chosen = candidates[order[0]]
    scores = [score_optimized(family.sites, chosen) for family in validation]
    tuning_rows = [
        {"stage": "coarse", **asdict(params), "mean_spearman": float(coarse_means[i])}
        for i, params in enumerate(coarse)
    ]
    tuning_rows.extend(
        {"stage": "refined_validation", **asdict(params), "mean_spearman": float(validation_means[i])}
        for i, params in enumerate(candidates)
    )
    chosen_segmentation = (float("-inf"), CONTEXT)
    for beta, minimum, maximum in product((2., 6., 10., 18.), (2, 3, 5), (None, 64)):
        context = replace(CONTEXT, beta=beta, min_segment=minimum, max_segment=maximum)
        f1 = []
        for family, raw in zip(validation, scores):
            if family.truth.boundaries:
                segments = contextual_scores(raw, family.sites.coverage, "segmentation", context).segments
                boundaries = tuple(segment.stop for segment in segments[:-1])
                f1.append(boundary_metrics(boundaries, family.truth.boundaries)[-1])
        objective = mean(f1)
        tuning_rows.append({"stage": "segmentation_validation", "beta": beta,
                            "min_segment": minimum, "max_segment": maximum, "boundary_f1": objective})
        if objective > chosen_segmentation[0] + 1e-12:
            chosen_segmentation = objective, context

    chosen_window = (float("-inf"), chosen_segmentation[1])
    for strength, radius in product((.05, .25, 1., 4.), (4, 8, 16)):
        context = replace(chosen_segmentation[1], smoothing=strength,
                          radii=tuple(r for r in (0, 1, 2, 4, 8, 16) if r <= radius))
        rho = []
        for family, raw in zip(validation, scores):
            if np.unique(family.truth.truth[family.sites.valid]).size > 1:
                smooth = contextual_scores(raw, family.sites.coverage, "scaw", context).scores
                valid = family.sites.valid
                correlation = float(correlations(smooth[valid], family.truth.truth[valid], rank=True))
                if np.isfinite(correlation):
                    rho.append(correlation)
        objective = mean(rho)
        tuning_rows.append({"stage": "window_validation", "smoothing": strength,
                            "radius": radius, "mean_spearman": objective})
        if objective > chosen_window[0] + 1e-12:
            chosen_window = objective, context

    regimes = [f.scenario for f in validation if
               np.unique(f.truth.truth[f.sites.valid]).size > 1]
    per_regime = {}
    for scenario in sorted(set(regimes)):
        means = _averages(validation_objectives[np.array(regimes) == scenario])
        best = _sorted(means, candidates)[0]
        per_regime[scenario] = {"best_parameters": asdict(candidates[best]),
                                "best_spearman": float(means[best]),
                                "global_selection_spearman": float(means[candidates.index(chosen)])}
    return {
        "reference": asdict(DEFAULT_PARAMETERS), "training_coarse_best": asdict(best_training),
        "rvne": asdict(chosen), "context": asdict(chosen_window[1]), "coarse_count": len(coarse),
        "training_grid_seconds": grid_seconds, "refined_count": len(candidates),
        "validation_spearman": float(validation_means[order[0]]),
        "reference_validation_spearman": float(validation_means[candidates.index(DEFAULT_PARAMETERS)]),
        "validation_boundary_f1": chosen_segmentation[0], "validation_context_spearman": chosen_window[0],
        "validation_by_scenario": per_regime,
    }, tuning_rows


def methods(family: Family, retuned: RVNEParameters, context: ContextParameters) -> dict:
    """Get one held-out score profile and optional segmentation per method."""
    sites = family.sites
    reference = score_reference(sites)
    raw = score_optimized(sites)
    tuned = score_optimized(sites, retuned)
    result = {"E0": (reference * sites.coverage, None), "E1": (raw * sites.coverage, None),
              "E7": (tuned * sites.coverage, None)}
    for label, method in CONTEXT_METHODS.items():
        fitted = contextual_scores(raw, sites.coverage, method)
        result[label] = fitted.scores, fitted
    fitted = contextual_scores(tuned, sites.coverage, "scaw", context)
    result["E8"] = fitted.scores, fitted
    result.update({label: (scores, None) for label, scores in
                   comparator_profiles(sites, family.truth.alignment).items()})
    return result


def calibration_thresholds(validation: list[Family], params: RVNEParameters,
                           context: ContextParameters) -> dict[str, float]:
    """Learn each binary decision threshold on validation scores only."""
    scores, labels = defaultdict(list), []
    for family in validation:
        valid = family.sites.valid
        labels.extend(family.truth.conserved[valid].tolist())
        for name, (prediction, _) in methods(family, params, context).items():
            scores[name].extend(np.round(prediction[valid], 12).tolist())
    return {name: select_threshold(values, labels) for name, values in scores.items()}


def evaluation_rows(inputs: list[Family], params: RVNEParameters,
                    context: ContextParameters, thresholds: dict[str, float]) -> list[dict]:
    """Summarize each independent test family before group aggregation."""
    rows = []
    for family in inputs:
        truth = family.truth
        for name, (scores, fitted) in methods(family, params, context).items():
            measures = site_metrics(scores, truth.truth, truth.conserved, thresholds[name])
            row = {"method": name, "seed": family.seed, "scenario": family.scenario,
                   "N": truth.alignment.size, "L": truth.alignment.length,
                   "settings": json.dumps(family.settings, sort_keys=True), **measures}
            valid = family.sites.valid
            if fitted is not None:
                widths = fitted.window_size[valid]
                row["median_window"] = float(np.median(widths)) if widths.size else float("nan")
                row["singleton_fraction"] = float(np.mean(widths == 1)) if widths.size else float("nan")
            if name in ("E5", "E6", "E8"):
                segments = fitted.segments
                predicted = tuple(s.stop for s in segments[:-1])
                lengths = [s.stop - s.start for s in segments]
                # Assign an entire segment on the coverage-penalized decision scale.
                regime = np.zeros(len(valid), dtype=bool)
                for segment in segments:
                    region = slice(segment.start, segment.stop)
                    mean_coverage = family.sites.coverage[region][valid[region]].mean()
                    regime[region] = round(segment.mean * mean_coverage, 12) >= thresholds[name]
                row.update(region_iou=region_iou(regime, truth.conserved, valid),
                           segments=len(segments), segment_lengths=",".join(map(str, lengths)))
                if truth.boundaries:
                    precision, recall, f1 = boundary_metrics(predicted, truth.boundaries)
                    row.update(boundary_precision=precision, boundary_recall=recall, boundary_f1=f1)
            rows.append(row)
    return rows


METRICS = ("scored_fraction", "spearman", "pearson", "auroc", "auprc", "sensitivity", "specificity",
           "boundary_precision", "boundary_recall", "boundary_f1", "region_iou", "segments",
           "median_window", "singleton_fraction")


def summarize(rows: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["method"], "all"].append(row)
        grouped[row["method"], row["scenario"]].append(row)
    summary = []
    for (method, scenario), group in sorted(grouped.items()):
        for metric in METRICS:
            values = [row[metric] for row in group if np.isfinite(row.get(metric, np.nan))]
            if values:
                summary.append(dict(method=method, scenario=scenario, metric=metric, n=len(values),
                                    mean=mean(values), sd=stdev(values) if len(values) > 1 else 0.0))
    return summary


def timed(function, repeats: int) -> dict:
    """Measure wall time separately from Python/NumPy traced allocation peak."""
    function()
    elapsed = []
    for _ in range(repeats):
        start = perf_counter()
        function()
        elapsed.append(perf_counter() - start)
    tracemalloc.start()
    try:
        function()
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    return {"median_seconds": median(elapsed), "sd_seconds": stdev(elapsed),
            "peak_traced_bytes": peak, "repeats": repeats}


def runtime_rows(repeats: int) -> list[dict]:
    rows = []
    for size, length in ((20, 100), (50, 100), (100, 100), (500, 100),
                         (50, 500), (50, 1000), (50, 5000)):
        alignment = simulate(size=size, length=length, seed=90100 + size + length).alignment
        sites = site_profiles(alignment)
        raw, optimized = score_reference(sites), score_optimized(sites)
        error = abs(raw - optimized)
        rows.append(dict(N=size, L=length, stage="equivalence", max_abs=float(error.max()),
                         mean_abs=float(error.mean()), spearman=float(correlations(raw, optimized, rank=True))))
        raw_reference = timed(lambda: score_reference(sites), repeats)
        fast = timed(lambda: score_optimized(sites), repeats)
        rows.extend([
            dict(N=size, L=length, stage="profiles", **timed(lambda: site_profiles(alignment), repeats)),
            dict(N=size, L=length, stage="reference", **raw_reference),
            dict(N=size, L=length, stage="optimized", **fast,
                 speedup=raw_reference["median_seconds"] / fast["median_seconds"]),
            dict(N=size, L=length, stage="context", **timed(
                lambda: contextual_scores(optimized, sites.coverage, "scaw"), repeats)),
        ])
    for length in (100, 500, 1000, 5000):
        sites = site_profiles(simulate(length=length, seed=90200 + length).alignment)
        grid = paper_parameter_grid()
        rows.append(dict(N=50, L=length, stage="grid_1331",
                         **timed(lambda: score_grid(sites, grid), repeats)))
    return rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    """Write rows without losing fields that are absent for some methods."""
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("tune", "benchmark", "all"), default="all")
    parser.add_argument("--parameters", type=Path, help="previously frozen tuning JSON")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--sweep", action="store_true", help="include held-out one-factor sweeps")
    args = parser.parse_args()
    if args.repeats < 2:
        parser.error("at least two timing repeats are required")
    if args.mode == "benchmark" and args.parameters is None:
        parser.error("benchmark mode requires --parameters")
    args.output.mkdir(parents=True, exist_ok=True)
    validation = families(BASE_SPECS, VALIDATION)
    if args.mode in ("tune", "all"):
        selected, tuning = select_parameters(families(BASE_SPECS, TRAINING), validation)
        write_tsv(args.output / "tuning.tsv", tuning)
    else:
        selected = json.loads(args.parameters.read_text(encoding="utf-8"))
    rvne = RVNEParameters(**selected["rvne"])
    context_args = selected["context"].copy()
    context_args["radii"] = tuple(context_args["radii"])
    context = ContextParameters(**context_args)
    thresholds = calibration_thresholds(validation, rvne, context)
    thresholds.update(selected.get("thresholds") or {})
    selected["thresholds"] = thresholds
    selected["training_seed_groups"] = TRAINING
    selected["validation_seed_groups"] = VALIDATION
    selected["test_seed_groups"] = HELD_OUT
    selected["python"] = platform.python_version()
    selected["numpy"] = np.__version__
    (args.output / "parameters.json").write_text(
        json.dumps(selected, allow_nan=False, indent=2) + "\n", encoding="utf-8")
    if args.mode != "tune":
        rows = evaluation_rows(families(BASE_SPECS, HELD_OUT), rvne, context, thresholds)
        write_tsv(args.output / "test_families.tsv", rows)
        write_tsv(args.output / "test_summary.tsv", summarize(rows))
        if args.sweep:
            rows = evaluation_rows(families(SWEEPS, (401, 402, 403)), rvne, context, thresholds)
            write_tsv(args.output / "sweep_families.tsv", rows)
            write_tsv(args.output / "sweep_summary.tsv", summarize(rows))
        write_tsv(args.output / "runtime.tsv", runtime_rows(args.repeats))
    print(f"Wrote {args.mode} outputs to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Evaluate current SCAW and its window and segmentation ablations."""

import argparse
from dataclasses import asdict
from functools import partial
import hashlib
import json
from pathlib import Path
import platform
import shutil
from time import perf_counter

import numpy as np

from data_generation.rvne_simulation import simulate, write_simulation
from data_generation.indelible import Options, find_executable
from msa_conservation.scoring.rvne import RVNEParameters, site_profiles
from msa_conservation.scoring.rvne_optimized import score_optimized
from msa_conservation.scoring.scaw import ContextParameters, contextual_scores, METHOD_VERSION
from .metrics import select_threshold
from .run_rvne_benchmark import BASE_SPECS, VALIDATION, HELD_OUT, Family, families, summarize, write_tsv
from .run_supplemental_benchmark import load_families
from .rvne_evaluation import boundary_metrics, region_iou, site_metrics


CODE = Path(__file__).resolve().parents[1]
METRICS = ("spearman", "auroc", "auprc", "sensitivity", "specificity", "scored_fraction",
           "boundary_f1", "boundary_f1_exact", "region_iou", "segments", "false_boundaries",
           "extra_boundaries", "missing_boundaries", "seconds")


def scoring_methods(raw, coverage):
    return {
        "v2": lambda: contextual_scores(raw, coverage, "scaw"),
        "v2_segment_mean": lambda: contextual_scores(raw, coverage, "segmentation"),
        "fixed_uniform_7": lambda: contextual_scores(raw, coverage, "uniform"),
    }


def calibrate(validation, rvne):
    scores, labels = {}, []
    sensitivity = []
    for family in validation:
        raw = score_optimized(family.sites, rvne)
        valid = family.sites.valid
        labels.extend(family.truth.conserved[valid].tolist())
        for name, score in scoring_methods(raw, family.sites.coverage).items():
            result = score()
            scores.setdefault(name, []).extend(np.round(result.scores[valid], 12).tolist())
        for scale in (.5, 1., 2.):
            result = contextual_scores(raw, family.sites.coverage, "scaw", ContextParameters(scale))
            boundaries = tuple(s.stop for s in result.segments[:-1])
            sensitivity.append({"method": str(scale), "scenario": family.scenario,
                                "seed": family.seed, "penalty_scale": scale,
                                **site_metrics(result.scores, family.truth.truth,
                                               family.truth.conserved, 0.),
                                "boundary_f1": (boundary_metrics(boundaries, family.truth.boundaries)[-1]
                                                if family.scenario != "gradual" else np.nan),
                                "false_boundaries": (len(boundaries) if not family.truth.boundaries
                                                     and family.scenario != "gradual" else np.nan)})
    return {name: select_threshold(values, labels) for name, values in scores.items()}, sensitivity


def evaluate(inputs, cohort, rvne, thresholds):
    rows, provenance = [], []
    for family in inputs:
        sites, truth = family.sites, family.truth
        raw = score_optimized(sites, rvne)
        provenance.append({"cohort": cohort, "seed": family.seed,
                           "alignment_sha256": hashlib.sha256(
                               "\n".join(truth.alignment.sequences).encode()).hexdigest(),
                           "rates_sha256": hashlib.sha256(truth.rates.tobytes()).hexdigest(),
                           "size": truth.alignment.size, "columns": truth.alignment.length})
        for name, function in scoring_methods(raw, sites.coverage).items():
            start = perf_counter()
            fitted = function()
            elapsed = perf_counter() - start
            row = {"cohort": cohort, "seed": family.seed, "method": name,
                   "scenario": family.scenario, "seconds": elapsed,
                   **site_metrics(fitted.scores, truth.truth, truth.conserved, thresholds[name])}
            if fitted.segments:
                predicted = tuple(s.stop for s in fitted.segments[:-1])
                precision, _, f1 = boundary_metrics(predicted, truth.boundaries)
                # A smooth truth has no discrete boundary targets; mark undefined.
                if family.scenario != "gradual":
                    row["boundary_f1"] = f1
                    row["boundary_f1_exact"] = boundary_metrics(predicted, truth.boundaries, 0)[-1]
                    row["false_boundaries"] = round(len(predicted) * (1 - precision))
                    row["extra_boundaries"] = max(0, len(predicted) - len(truth.boundaries))
                    row["missing_boundaries"] = max(0, len(truth.boundaries) - len(predicted))
                regime = np.zeros(len(raw), dtype=bool)
                for segment in fitted.segments:
                    region = slice(segment.start, segment.stop)
                    mean_coverage = sites.coverage[region][sites.valid[region]].mean()
                    regime[region] = round(segment.mean * mean_coverage, 12) >= thresholds[name]
                row["region_iou"] = region_iou(regime, truth.conserved, sites.valid)
                row["segments"] = len(fitted.segments)
            rows.append(row)
    return rows, provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new directory")
    parser.add_argument("--fresh-indelible", action="store_true", help="include the 12 dedicated INDELible cases")
    parser.add_argument("--manifest", type=Path, action="append", default=[],
                        help="include a generated supplemental suite manifest (repeatable)")
    args = parser.parse_args()
    try:
        executable = find_executable()
        supplemental = [family for path in args.manifest for family in load_families(path)]
        if len({family.seed for family in supplemental}) != len(supplemental):
            raise ValueError("supplemental manifests repeat a family seed")
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    args.output.mkdir(parents=True, exist_ok=False)
    # Freeze the RVNE values used for this evaluation so removing old SCAW
    # artifacts does not change the profiles or the benchmark split.
    rvne = RVNEParameters(r=0., s1=.1, s2=0.)
    generate = partial(simulate, executable=executable)
    validation = families(BASE_SPECS, VALIDATION, generate)
    thresholds, sensitivity = calibrate(validation, rvne)
    # c=1 and radius=3 were specified before test evaluation. Sensitivity is
    # diagnostic only, on validation families; it does not select the default.
    parameters = {"context_version": METHOD_VERSION, "rvne": asdict(rvne),
                  "context": asdict(ContextParameters()), "thresholds": thresholds,
                  "default_selection": "c=1, radius=3, center=0.5 fixed before test evaluation",
                  "rvne_source": "Fixed RVNE values from the original training/validation split",
                  "threshold_source": "24 INDELible 1.03 validation families; calibrated per method",
                  "implementation_sha256": hashlib.sha256(
                      (CODE / "rvne/src/msa_conservation/scoring/scaw.py").read_bytes()).hexdigest(),
                  "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  "python": platform.python_version(), "numpy": np.__version__}
    (args.output / "parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    write_tsv(args.output / "validation_sensitivity.tsv", sensitivity)
    write_tsv(args.output / "validation_summary.tsv", summarize(sensitivity, METRICS))
    cohorts = {"indelible_core": families(BASE_SPECS, HELD_OUT, generate)}
    if supplemental:
        cohorts["indelible_supplemental"] = supplemental
    if args.fresh_indelible:
        fresh = []
        for seed, scenario, motif, indel in (
            (81101, "conserved", 3, 0.), (81102, "conserved-variable-blocks", 3, 0.),
            (81103, "single-motif", 3, 0.), (81104, "single-motif", 5, 0.),
            (81105, "three-rate-blocks", 3, 0.), (81106, "three-rate-blocks", 3, .03),
            (81201, "conserved", 3, 0.), (81202, "conserved-variable-blocks", 3, 0.),
            (81203, "single-motif", 3, 0.), (81204, "single-motif", 5, 0.),
            (81205, "three-rate-blocks", 3, 0.), (81206, "three-rate-blocks", 3, .03),
        ):
            folder = args.output / "fresh_inputs" / str(seed)
            simulation = generate(seed=seed, size=20, length=120, scenario=scenario,
                                  motif_length=motif, work_dir=folder,
                                  options=Options(insertion_rate=indel, deletion_rate=indel))
            write_simulation(simulation, folder / "alignment.fasta", folder / "truth.tsv", folder)
            fresh.append(Family(seed, scenario, {"motif_length": motif, "indel_rate": indel},
                                simulation, site_profiles(simulation.alignment)))
        cohorts["fresh_indelible"] = fresh
    all_rows, provenance, summaries = [], [], []
    for cohort, inputs in cohorts.items():
        rows, origin = evaluate(inputs, cohort, rvne, thresholds)
        all_rows.extend(rows)
        provenance.extend(origin)
        summaries.extend({"cohort": cohort, **row} for row in summarize(rows, METRICS))
    write_tsv(args.output / "families.tsv", all_rows)
    write_tsv(args.output / "summary.tsv", summaries)
    write_tsv(args.output / "inputs.tsv", provenance)
    # Context-only timing, a shared noisy step profile, three measured repetitions.
    runtime = []
    rng = np.random.default_rng(912)
    for length in (240, 1000, 5000):
        raw = np.repeat([.2, .7, .3, 1.], (length + 3) // 4)[:length]
        raw = raw + rng.normal(0, .03, length)
        coverage = np.ones(length)
        function = scoring_methods(raw, coverage)["v2"]
        function()
        elapsed = []
        for _ in range(3):
            start = perf_counter()
            function()
            elapsed.append(perf_counter() - start)
        runtime.append({"method": "v2", "columns": length, "repeats": 3,
                        "median_seconds": float(np.median(elapsed))})
    write_tsv(args.output / "runtime.tsv", runtime)
    if args.fresh_indelible:
        shutil.make_archive(str(args.output / "fresh_inputs"), "zip", args.output / "fresh_inputs")
    print(json.dumps({"families": {key: len(value) for key, value in cohorts.items()},
                      "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()

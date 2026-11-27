"""Build analysis tables and exact plotting data from the frozen benchmark."""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from benchmarks.run_rvne_benchmark import methods
from benchmarks.run_supplemental_benchmark import load_families
from benchmarks.rvne_evaluation import average_precision
from msa_conservation.scoring.rvne import RVNEParameters
from msa_conservation.scoring.scaw import ContextParameters


OUTPUT = Path(__file__).resolve().parent
RVNE = OUTPUT.parent / "rvne"
DATASETS = OUTPUT.parent / "data_generation" / "datasets"
TABLES = OUTPUT / "tables"
SUPPLEMENTAL = RVNE / "results" / "supplemental_run" / "families.tsv"
REFERENCE = RVNE / "results" / "reference_run"
METHODS = [f"E{i}" for i in range(9)] + [
    "Shannon", "KL", "JSD", "BLOSUM62", "ClustalW consensus", "AL2CO variance",
]
METRICS = {
    "A": ("variable_fpr", "false_boundaries", "scored_fraction"),
    "B": ("spearman", "auroc", "auprc", "variable_recall", "variable_specificity",
          "variable_site_iou", "boundary_f1"),
    "C": ("spearman", "auroc", "auprc", "sensitivity", "specificity", "boundary_f1"),
    "D": ("spearman", "auroc", "auprc", "variable_recall", "variable_specificity",
          "variable_site_iou", "region_iou", "boundary_f1", "medium_as_variable",
          "high_as_variable", "auroc_low_vs_medium", "auroc_medium_vs_high"),
}
SETTING = {"A": "divergence", "B": "variable_fraction", "C": "motif_length", "D": "variable_fraction"}


def read(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(name, rows):
    with (TABLES / f"{name}.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def numeric(row, metric):
    return float(row.get(metric) or "nan")


def summary_rows(groups):
    out = []
    for keys, values in groups.items():
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        n = len(finite)
        out.append(dict((*keys, ("n", n), ("mean", finite.mean() if n else ""),
                         ("sd", finite.std(ddof=1) if n > 1 else ""))))
    return out


def supplemental_stats(rows):
    by_pattern, by_setting = defaultdict(list), defaultdict(list)
    for row in rows:
        pattern, method = row["scenario"], row["method"]
        setting = SETTING[pattern]
        level = json.loads(row["settings"])[setting]
        for metric in METRICS[pattern]:
            value = (1 - numeric(row, "variable_specificity") if metric == "variable_fpr"
                     else numeric(row, metric))
            keys = (("pattern", pattern), ("method", method), ("metric", metric))
            by_pattern[keys].append(value)
            by_setting[keys + (("setting", setting), ("level", level))].append(value)
    write("table_supplemental", summary_rows(by_pattern))
    write("table_by_setting", summary_rows(by_setting))


def paired_stats(rows):
    indexed = {(int(row["seed"]), row["method"]): row for row in rows if row["scenario"] == "D"}
    seeds = sorted({seed for seed, _ in indexed})
    strata = defaultdict(list)
    for index, seed in enumerate(seeds):
        strata[json.loads(indexed[seed, "E8"]["settings"])["variable_fraction"]].append(index)
    generator = np.random.default_rng(41027)
    samples = np.concatenate([
        generator.choice(indices, size=(10000, len(indices)), replace=True) for indices in strata.values()
    ], axis=1)
    effects = []
    for first, second in (("E8", "E7"), ("E8", "E5"), ("E6", "E4")):
        for metric in ("variable_recall", "medium_as_variable", "variable_site_iou", "spearman"):
            delta = np.array([numeric(indexed[s, first], metric) - numeric(indexed[s, second], metric)
                              for s in seeds])
            low, high = np.quantile(delta[samples].mean(axis=1), [.025, .975])
            wins = int((delta > 1e-12).sum())
            ties = int((np.abs(delta) <= 1e-12).sum())
            effects.append(dict(first=first, second=second, metric=metric, n=len(delta),
                                mean_delta=delta.mean(), ci_low=low, ci_high=high,
                                wins=wins, ties=ties, losses=len(delta)-wins-ties))
    write("table_d_paired", effects)


def reference_stats():
    core, gap = defaultdict(list), defaultdict(list)
    for row in read(REFERENCE / "test_families.tsv"):
        for metric in ("spearman", "auroc", "auprc", "sensitivity", "specificity", "boundary_f1", "region_iou"):
            core[(("method", row["method"]), ("metric", metric))].append(numeric(row, metric))
    for row in read(REFERENCE / "sweep_families.tsv"):
        level = json.loads(row["settings"]).get("gaps")
        if level is None:
            continue
        for metric in ("scored_fraction", "spearman", "auroc", "boundary_f1", "region_iou"):
            keys = (("method", row["method"]), ("nominal_gap", level), ("metric", metric))
            gap[keys].append(numeric(row, metric))
    write("table_reference_core", summary_rows(core))
    write("table_reference_gap", summary_rows(gap))


def curve_rows(method, scores, labels):
    order = np.argsort(-scores, kind="stable")
    scores, labels = scores[order], labels[order]
    endpoints = np.r_[np.flatnonzero(scores[1:] != scores[:-1]), len(scores)-1]
    tp = np.r_[0, np.cumsum(labels)[endpoints]]
    fp = np.r_[0, endpoints+1] - tp
    tpr, fpr = tp / tp[-1], fp / fp[-1]
    precision = np.divide(tp, tp+fp, out=np.ones(len(tp)), where=tp+fp > 0)
    auroc = np.trapezoid(tpr, fpr)
    ap = average_precision(scores, labels)
    assert np.isclose(ap, np.sum(np.diff(tpr) * precision[1:]))
    return [dict(method=method, fpr=x, recall=r, precision=p, auroc=auroc, ap=ap,
                 prevalence=np.mean(labels)) for x, r, p in zip(fpr, tpr, precision)]


def site_level(rows):
    config = json.loads((REFERENCE / "parameters.json").read_text(encoding="utf-8"))
    context = ContextParameters(**{**config["context"], "radii": tuple(config["context"]["radii"])})
    parameters = RVNEParameters(**config["rvne"])
    thresholds = config["thresholds"]
    manifests = [DATASETS / pattern / "manifest.tsv" for pattern in ("patterns_abc", "patterns_d")]
    families = {family.seed: family for manifest in manifests for family in load_families(manifest)}
    selections = {50124: ("E2", "E3", "E8"), 50140: ("E5", "E7", "E8")}
    curve_methods = ("E0", "E5", "E7", "E8", "Shannon", "AL2CO variance")
    curves, profile = [], []
    d_scores = {method: [] for method in curve_methods}
    d_labels = {method: [] for method in curve_methods}
    lookup = {(int(row["seed"]), row["method"]): row for row in rows}
    for seed, family in families.items():
        if seed not in selections and family.scenario != "D":
            continue
        profiles = methods(family, parameters, context)
        for method in selections.get(seed, ()):
            scores = profiles[method][0]
            for i, (score, rate) in enumerate(zip(scores, family.truth.rates), 1):
                profile.append(dict(seed=seed, pattern=family.scenario, method=method, position=i,
                                    rate=rate, score=score, threshold=thresholds[method],
                                    margin=np.round(score, 12)-thresholds[method]))
        if family.scenario == "D":
            for method in curve_methods:
                scores = np.round(profiles[method][0], 12)
                valid = np.isfinite(scores)
                labels = ~family.truth.conserved[valid]
                values = -scores[valid]
                assert np.isclose(np.mean(values[labels] > -thresholds[method]),
                                  numeric(lookup[seed, method], "variable_recall"))
                d_scores[method].extend(values)
                d_labels[method].extend(labels)
    for method in curve_methods:
        curves.extend(curve_rows(method, np.array(d_scores[method]), np.array(d_labels[method])))
    write("figure_position_profiles", profile)
    write("figure_d_curves", curves)


def main():
    TABLES.mkdir(exist_ok=True)
    rows = read(SUPPLEMENTAL)
    expected = {pattern: 12 * len(METHODS) for pattern in SETTING}
    if Counter(row["scenario"] for row in rows) != expected:
        raise ValueError("expected 12 supplemental families per pattern and 15 methods per family")
    supplemental_stats(rows)
    paired_stats(rows)
    reference_stats()
    site_level(rows)


if __name__ == "__main__":
    main()

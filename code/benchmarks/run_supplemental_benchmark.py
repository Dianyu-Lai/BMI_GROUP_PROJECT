"""Evaluate aligned FASTA and matched truth TSVs with frozen RVNE settings."""

import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path

import numpy as np

from data_generation.rvne_simulation import Simulation
from msa_conservation.msa import read_fasta
from msa_conservation.scoring.rvne import RVNEParameters, site_profiles
from msa_conservation.scoring.scaw import ContextParameters
from .metrics import roc_auc
from .run_rvne_benchmark import METRICS, Family, evaluation_rows, summarize, write_tsv
from .rvne_evaluation import region_iou


DIAGNOSTICS = (
    "variable_recall", "variable_specificity", "variable_site_iou", "false_boundaries",
    "auroc_low_vs_medium", "auroc_medium_vs_high", "medium_as_variable", "high_as_variable",
)
PATTERN_SETTINGS = {
    "A": "divergence", "B": "variable_fraction", "C": "motif_length", "D": "variable_fraction",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_families(manifest: Path) -> list[Family]:
    """Load the supplied alignments and translate one-based boundary markers."""
    families = []
    for row in read_tsv(manifest):
        fasta = manifest.parent / row["fasta"]
        truth_path = manifest.parent / row["truth"]
        with fasta.open(encoding="utf-8") as handle:
            alignment = read_fasta(handle)
        if (alignment.size, alignment.length) != (int(row["size"]), int(row["length"])):
            raise ValueError(f"alignment dimensions disagree with {manifest}: {fasta}")
        truth = read_tsv(truth_path)
        positions = [int(site["position"]) for site in truth]
        if positions != list(range(1, alignment.length + 1)):
            raise ValueError(f"truth must cover every aligned position in order: {truth_path}")
        rates = np.array([float(site["rate"]) for site in truth])
        conservation = np.array([float(site["conservation"]) for site in truth])
        labels = np.array([int(site["conserved"]) for site in truth])
        markers = np.array([int(site["boundary_before"]) for site in truth])
        boundaries = tuple(np.flatnonzero(markers).tolist())
        if (not np.all(np.isfinite(rates) & (rates >= 0)) or
                not np.all((0 <= conservation) & (conservation <= 1)) or
                not np.all(np.isin(labels, (0, 1))) or
                not np.all(np.isin(markers, (0, 1))) or
                boundaries != tuple((np.flatnonzero(np.diff(rates)) + 1).tolist())):
            raise ValueError(f"invalid rates, conservation values, labels or boundaries: {truth_path}")
        simulation = Simulation(alignment, rates, conservation, labels.astype(bool), boundaries)
        families.append(Family(int(row["seed"]), row["pattern"], row, simulation, site_profiles(alignment)))
    if not families:
        raise ValueError(f"manifest contains no alignments: {manifest}")
    return families


def compare_rates(scores: np.ndarray, rates: np.ndarray, lower: float, higher: float) -> float:
    selected = (rates == lower) | (rates == higher)
    labels = rates[selected] == lower
    if not labels.any() or labels.all():
        return float("nan")
    return roc_auc(scores[selected].tolist(), labels.tolist())


def pattern_diagnostics(family: Family, scores: np.ndarray, threshold: float, row: dict) -> dict:
    """Report variable-site metrics and distinguish the adjacent rates in D."""
    valid = np.isfinite(scores)
    rounded = np.round(scores, 12)
    variable = rounded < threshold
    extra = {
        "variable_recall": row["specificity"],
        "variable_specificity": row["sensitivity"],
        "variable_site_iou": region_iou(variable, ~family.truth.conserved, valid),
    }
    if not family.truth.boundaries and "segments" in row:
        extra["false_boundaries"] = max(0, row["segments"] - 1)
    if family.scenario == "D":
        low, medium, high = np.unique(family.truth.rates)
        scores_valid, rates_valid = rounded[valid], family.truth.rates[valid]
        extra["auroc_low_vs_medium"] = compare_rates(scores_valid, rates_valid, low, medium)
        extra["auroc_medium_vs_high"] = compare_rates(scores_valid, rates_valid, medium, high)
        for metric, rate in (("medium_as_variable", medium), ("high_as_variable", high)):
            selected = valid & (family.truth.rates == rate)
            extra[metric] = float(np.mean(variable[selected])) if selected.any() else float("nan")
    return extra


def setting_summaries(rows: list[dict]) -> list[dict]:
    settings = []
    for row in rows:
        entry = json.loads(row["settings"])
        name = PATTERN_SETTINGS[row["scenario"]]
        settings.append({**row, "scenario": f"{row['scenario']} {name}={entry[name]}"})
    return [row for row in summarize(settings, METRICS + DIAGNOSTICS) if row["scenario"] != "all"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, action="append", required=True)
    parser.add_argument("--parameters", type=Path, required=True, help="frozen parameters and thresholds")
    parser.add_argument("--output", type=Path, required=True, help="new output directory")
    args = parser.parse_args()
    try:
        selected = json.loads(args.parameters.read_text(encoding="utf-8"))
        rvne = RVNEParameters(**selected["rvne"])
        context = ContextParameters(**{**selected["context"], "radii": tuple(selected["context"]["radii"])})
        thresholds = selected["thresholds"]
        families = [family for manifest in args.manifest for family in load_families(manifest)]
        if len({family.seed for family in families}) != len(families):
            raise ValueError("input manifests repeat a family seed")
        rows = evaluation_rows(families, rvne, context, thresholds, diagnostics=pattern_diagnostics)
        summary = summarize(rows, METRICS + DIAGNOSTICS)
        by_setting = setting_summaries(rows)
        args.output.mkdir(parents=True)
        write_tsv(args.output / "families.tsv", rows)
        write_tsv(args.output / "summary.tsv", summary)
        write_tsv(args.output / "settings_summary.tsv", by_setting)
        metadata = {
            "manifests": [path.as_posix() for path in args.manifest],
            "parameters": args.parameters.as_posix(),
            "parameters_sha256": hashlib.sha256(args.parameters.read_bytes()).hexdigest(),
            "family_count": len(families),
            "method_count": len({row["method"] for row in rows}),
            "python": platform.python_version(),
            "numpy": np.__version__,
        }
        (args.output / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(f"Scored {metadata['family_count']} families with {metadata['method_count']} methods in {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

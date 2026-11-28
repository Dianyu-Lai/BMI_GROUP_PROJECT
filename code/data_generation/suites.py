"""Generate supplemental A/B/C and three-rate D data with optional rate-specific gaps."""

from collections import defaultdict
import csv
import json
from pathlib import Path
from statistics import mean

from .indelible import Options
from .rvne_simulation import DEFAULT_DIVERGENCE, simulate, write_simulation


D_PATTERN = (("D", "three-rate-blocks", "variable_fraction", (0.25, 0.4, 0.6)),)
D_GAP_TARGETS = (0.03, 0.12, 0.25)
SUITES = {
    "abc": (50100, (
        ("A", "conserved", "divergence", (0.35, 0.7, 1.4)),
        ("B", "conserved-variable-blocks", "variable_fraction", (0.25, 0.4, 0.6)),
        ("C", "single-motif", "motif_length", (3, 5, 9)),
    )),
    "d": (50136, D_PATTERN),
    "d-gaps": (50148, D_PATTERN),
}


def generate_suite(
    output: Path, suite: str, *, size: int, length: int, replicates: int,
    start_seed: int | None, gap_targets: tuple[float, ...], executable: Path,
    options: Options, boundary_width: int,
) -> int:
    """Write a benchmark manifest and its INDELible families to a new directory."""
    default_seed, patterns = SUITES[suite]
    seed = default_seed if start_seed is None else start_seed
    with_rate_gaps = suite == "d-gaps"
    columns = ("pattern", "scenario", "seed", "size", "length", "ancestral_length", "divergence",
               "variable_fraction", "motif_length", "fasta", "truth", "metadata", "binary_metrics_reason")
    if with_rate_gaps:
        columns += ("reference_seed", "target_gap_low", "target_gap_medium", "target_gap_high")
    if replicates < 1 or size < 2 or length < 12 or boundary_width < 0:
        raise ValueError("require positive replicates, at least two taxa, 12 ancestral sites and nonnegative boundary width")
    if with_rate_gaps and not all(0 <= gap < 1 for gap in gap_targets):
        raise ValueError("all gap targets must be in [0, 1)")
    output.mkdir(parents=True)
    rows = []
    distances = defaultdict(lambda: defaultdict(list))
    for pattern, scenario, parameter, values in patterns:
        for value in values:
            settings = {"size": size, "length": length, "divergence": DEFAULT_DIVERGENCE}
            settings[parameter] = value
            for _ in range(replicates):
                stem = f"{pattern.lower()}_{seed}"
                run_path = Path("runs") / stem
                fasta = run_path / "alignment.fasta"
                truth = run_path / "ground_truth.tsv"
                simulation_seed = seed
                targets = None
                if with_rate_gaps:
                    simulation_seed = seed - (SUITES["d-gaps"][0] - SUITES["d"][0])
                    targets = tuple(gap_targets)
                simulation = simulate(**settings, scenario=scenario, seed=simulation_seed,
                                      gap_by_rate=targets, executable=executable, options=options,
                                      work_dir=output / run_path, boundary_width=boundary_width)
                write_simulation(simulation, output / fasta, output / truth, output / run_path)
                reason = simulation.metadata["binary_metrics_reason"]
                row = (pattern, scenario, seed, size, simulation.alignment.length, length,
                       settings["divergence"], settings.get("variable_fraction", ""),
                       settings.get("motif_length", ""), fasta.as_posix(), truth.as_posix(),
                       (run_path / "metadata.json").as_posix(),
                       reason or "")
                if with_rate_gaps:
                    row += (simulation_seed, *targets)
                rows.append(row)
                for part in simulation.metadata["partitions"]:
                    p_distance = part["statistics"]["mean_pairwise_p_distance"]
                    if p_distance is not None:
                        distances[pattern][part["rate"]].append(p_distance)
                seed += 1
    with (output / "manifest.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
    variability = {}
    for pattern, by_rate in distances.items():
        rates = {rate: mean(by_rate[rate]) for rate in sorted(by_rate)}
        values = list(rates.values())
        variability[pattern] = {"mean_partition_p_distance_by_rate": rates,
                                "ordered": all(left < right for left, right in zip(values, values[1:]))
                                if len(values) > 1 else None}
    (output / "validation.json").write_text(json.dumps(variability, indent=2) + "\n", encoding="utf-8")
    if any(row[0] == "A" for row in rows):
        print("Pattern A has one designed class. AUROC and AUPRC are undefined for its families.")
    return len(rows)

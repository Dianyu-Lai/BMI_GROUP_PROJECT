"""Generate supplemental A/B/C and three-rate D data with optional rate-specific gaps."""

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
from statistics import mean

from .indelible import find_executable
from .rvne_simulation import add_engine_arguments, engine_options, simulate, write_simulation


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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_engine_arguments(parser)
    parser.add_argument("--output", type=Path, required=True, help="new dataset directory")
    parser.add_argument("--suite", choices=SUITES, default="abc")
    parser.add_argument("--size", type=int, default=50)
    parser.add_argument("--length", type=int, default=240, help="ancestral sequence length")
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument("--start-seed", type=int, help="override the suite's default first family seed")
    parser.add_argument("--gap-targets", type=float, nargs=3, default=D_GAP_TARGETS,
                        metavar=("LOW", "MEDIUM", "HIGH"))
    args = parser.parse_args()
    default_seed, patterns = SUITES[args.suite]
    seed = default_seed if args.start_seed is None else args.start_seed
    with_rate_gaps = args.suite == "d-gaps"
    columns = ("pattern", "scenario", "seed", "size", "length", "ancestral_length", "divergence",
               "variable_fraction", "motif_length", "fasta", "truth", "metadata", "binary_metrics_reason")
    if with_rate_gaps:
        columns += ("reference_seed", "target_gap_low", "target_gap_medium", "target_gap_high")
    try:
        if args.replicates < 1 or args.size < 2 or args.length < 12 or args.boundary_width < 0:
            raise ValueError("require positive replicates, at least two taxa, 12 ancestral sites and nonnegative boundary width")
        if with_rate_gaps and not all(0 <= gap < 1 for gap in args.gap_targets):
            raise ValueError("all gap targets must be in [0, 1)")
        executable = find_executable(args.indelible_bin)
        options = engine_options(args)
        args.output.mkdir(parents=True)
        rows = []
        distances = defaultdict(lambda: defaultdict(list))
        for pattern, scenario, parameter, values in patterns:
            for value in values:
                settings = {"size": args.size, "length": args.length, "divergence": 0.7}
                settings[parameter] = value
                for _ in range(args.replicates):
                    stem = f"{pattern.lower()}_{seed}"
                    fasta, truth = f"{stem}.fasta", f"{stem}.truth.tsv"
                    run_path = Path("runs") / stem
                    simulation_seed = seed
                    targets = None
                    if with_rate_gaps:
                        simulation_seed = seed - (SUITES["d-gaps"][0] - SUITES["d"][0])
                        targets = tuple(args.gap_targets)
                    simulation = simulate(**settings, scenario=scenario, seed=simulation_seed,
                                          gap_by_rate=targets, executable=executable, options=options,
                                          work_dir=args.output / run_path, boundary_width=args.boundary_width)
                    write_simulation(simulation, args.output / fasta, args.output / truth, args.output / run_path)
                    reason = simulation.metadata["binary_metrics_reason"]
                    row = (pattern, scenario, seed, args.size, simulation.alignment.length, args.length,
                           settings["divergence"], settings.get("variable_fraction", ""),
                           settings.get("motif_length", ""), fasta, truth, (run_path / "metadata.json").as_posix(),
                           reason or "")
                    if with_rate_gaps:
                        row += (simulation_seed, *targets)
                    rows.append(row)
                    for part in simulation.metadata["partitions"]:
                        p_distance = part["statistics"]["mean_pairwise_p_distance"]
                        if p_distance is not None:
                            distances[pattern][part["rate"]].append(p_distance)
                    seed += 1
        with (args.output / "manifest.tsv").open("x", encoding="utf-8", newline="") as handle:
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
        (args.output / "validation.json").write_text(json.dumps(variability, indent=2) + "\n", encoding="utf-8")
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    if any(row[0] == "A" for row in rows):
        print("Pattern A has one designed class. AUROC and AUPRC are undefined for its families.")
    print(f"Wrote {len(rows)} families to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

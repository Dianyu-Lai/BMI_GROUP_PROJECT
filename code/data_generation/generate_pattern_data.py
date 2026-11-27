"""Generate supplemental A/B/C and three-rate D data with optional rate-specific gaps."""

import argparse
import csv
from pathlib import Path

from .rvne_simulation import simulate, write_simulation


D_PATTERN = (("D", "three-rate-blocks", "variable_fraction", (0.25, 0.4, 0.6)),)
D_GAP_PROBABILITIES = (0.03, 0.12, 0.25)
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
    parser.add_argument("--output", type=Path, required=True, help="new dataset directory")
    parser.add_argument("--suite", choices=SUITES, default="abc")
    args = parser.parse_args()
    seed, patterns = SUITES[args.suite]
    with_rate_gaps = args.suite == "d-gaps"
    columns = ("pattern", "scenario", "seed", "size", "length", "divergence",
               "variable_fraction", "motif_length", "fasta", "truth")
    if with_rate_gaps:
        columns += ("substitution_seed", "gap_seed", "gap_low", "gap_medium", "gap_high")
    try:
        args.output.mkdir(parents=True)
        with (args.output / "manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            for pattern, scenario, parameter, values in patterns:
                for value in values:
                    settings = {"size": 50, "length": 240, "divergence": 0.7}
                    settings[parameter] = value
                    for _ in range(4):
                        stem = f"{pattern.lower()}_{seed}"
                        fasta, truth = f"{stem}.fasta", f"{stem}.truth.tsv"
                        simulation_options = {**settings, "scenario": scenario, "seed": seed}
                        row = (pattern, scenario, seed, settings["size"], settings["length"],
                               settings["divergence"], settings.get("variable_fraction", ""),
                               settings.get("motif_length", ""), fasta, truth)
                        if with_rate_gaps:
                            substitution_seed = SUITES["d"][0] + seed - SUITES["d-gaps"][0]
                            simulation_options.update(seed=substitution_seed, gap_seed=seed,
                                                      gap_by_rate=D_GAP_PROBABILITIES)
                            row += (substitution_seed, seed, *D_GAP_PROBABILITIES)
                        simulation = simulate(**simulation_options)
                        write_simulation(simulation, args.output / fasta, args.output / truth)
                        writer.writerow(row)
                        seed += 1
    except OSError as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

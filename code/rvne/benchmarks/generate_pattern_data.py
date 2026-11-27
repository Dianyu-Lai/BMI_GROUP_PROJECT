"""Generate the supplemental A/B/C or three-rate D conservation-pattern data."""

import argparse
import csv
from pathlib import Path

from .rvne_simulation import simulate, write_simulation


SUITES = {
    "abc": (50100, (
        ("A", "conserved", "divergence", (0.35, 0.7, 1.4)),
        ("B", "conserved-variable-blocks", "variable_fraction", (0.25, 0.4, 0.6)),
        ("C", "single-motif", "motif_length", (3, 5, 9)),
    )),
    "d": (50136, (
        ("D", "three-rate-blocks", "variable_fraction", (0.25, 0.4, 0.6)),
    )),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new dataset directory")
    parser.add_argument("--suite", choices=SUITES, default="abc")
    args = parser.parse_args()
    seed, patterns = SUITES[args.suite]
    try:
        args.output.mkdir(parents=True)
        with (args.output / "manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("pattern", "scenario", "seed", "size", "length", "divergence",
                             "variable_fraction", "motif_length", "fasta", "truth"))
            for pattern, scenario, parameter, values in patterns:
                for value in values:
                    settings = {"size": 50, "length": 240, "divergence": 0.7}
                    settings[parameter] = value
                    for _ in range(4):
                        stem = f"{pattern.lower()}_{seed}"
                        fasta, truth = f"{stem}.fasta", f"{stem}.truth.tsv"
                        simulation = simulate(scenario=scenario, seed=seed, **settings)
                        write_simulation(simulation, args.output / fasta, args.output / truth)
                        writer.writerow((
                            pattern, scenario, seed, settings["size"], settings["length"],
                            settings["divergence"], settings.get("variable_fraction", ""),
                            settings.get("motif_length", ""), fasta, truth,
                        ))
                        seed += 1
    except OSError as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

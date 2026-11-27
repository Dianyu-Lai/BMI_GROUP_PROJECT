"""Score an aligned protein FASTA file and write tab-separated site scores."""

import argparse
import csv
import sys
from dataclasses import astuple, fields
from pathlib import Path
from typing import Sequence

from .msa import henikoff_weights, read_fasta
from .scoring.fusion import score_alignment
from .scoring.physicochemical import calibrated_blosum62


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("alignment", type=Path, help="aligned protein FASTA file")
    parser.add_argument("--weights", choices=("uniform", "henikoff"), default="uniform")
    parser.add_argument("--gamma", type=float, default=1.0, help="fusion coefficient (default: 1)")
    parser.add_argument("--gap-penalty", type=float, default=0.0, help="coefficient for missing nongap coverage")
    parser.add_argument("--aaindex-eta", type=float, default=0.0, help="AAindex tilt strength (default: 0)")
    parser.add_argument(
        "--reference-baselines", action="store_true", help="also score KL, JSD, BLOSUM62 SOP and adjusted entropy"
    )
    args = parser.parse_args(argv)

    try:
        with args.alignment.open(encoding="utf-8-sig") as file:
            alignment = read_fasta(file)
        weights = henikoff_weights(alignment) if args.weights == "henikoff" else None
        scores = score_alignment(
            alignment,
            weights=weights,
            gamma=args.gamma,
            gap_penalty=args.gap_penalty,
            matrix=calibrated_blosum62(args.aaindex_eta),
            reference_baselines=args.reference_baselines,
        )
        first = next(scores)
        writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        writer.writerow([field.name for field in fields(first)])
        writer.writerow(astuple(first))
        writer.writerows(astuple(site) for site in scores)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

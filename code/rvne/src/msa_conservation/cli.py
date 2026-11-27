"""Score aligned proteins with RVNE or run the real-protein alignment workflow."""

import argparse
from typing import Sequence


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("score", "real-data"))
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    if args.command == "score":
        from .rvne_cli import main as score_main

        return score_main(args.arguments)

    from .real_data import main as real_data_main

    return real_data_main(args.arguments)


if __name__ == "__main__":
    raise SystemExit(main())

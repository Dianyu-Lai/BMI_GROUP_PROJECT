"""Command line for single-family and suite generation."""

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

from .indelible import LengthModel, Options, find_executable, run_partition
from .rvne_simulation import (
    DEFAULT_DIVERGENCE, DEFAULT_LENGTH, DEFAULT_SIZE, Region, SCENARIOS, SUPPLEMENTAL_SCENARIOS,
    simulate, write_simulation,
)
from .suites import D_GAP_TARGETS, SUITES, generate_suite


def _add_engine_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--indelible-bin", type=Path)
    defaults = Options()
    length_default = (defaults.insertion_length.name, *map(str, defaults.insertion_length.parameters))
    parser.add_argument("--model", choices=("WAG", "LG", "JTT"), default=defaults.model)
    parser.add_argument("--method", type=int, choices=(1, 2), default=defaults.method)
    parser.add_argument("--insertion-rate", type=float, default=defaults.insertion_rate)
    parser.add_argument("--deletion-rate", type=float, default=defaults.deletion_rate)
    parser.add_argument("--insertion-length", nargs="+", default=length_default, metavar="MODEL_OR_PARAMETER")
    deletion = defaults.deletion_length
    parser.add_argument("--deletion-length", nargs="+", default=(deletion.name, *map(str, deletion.parameters)),
                        metavar="MODEL_OR_PARAMETER")
    parser.add_argument("--timeout", type=float, default=defaults.timeout, help="seconds per simulator process")
    parser.add_argument("--boundary-width", type=int, default=0, help="mask this many columns on each side")


def _engine_options(args: argparse.Namespace) -> Options:
    return Options(args.model, args.method, args.insertion_rate, args.deletion_rate,
                   LengthModel(args.insertion_length[0], tuple(map(float, args.insertion_length[1:]))),
                   LengthModel(args.deletion_length[0], tuple(map(float, args.deletion_length[1:]))), args.timeout)


def _read_regions(path: Path, options: Options) -> tuple[Region, ...]:
    entries = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(entries, list) or not entries:
        raise ValueError("region config must be a nonempty JSON array")
    regions = []
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("each region must be a JSON object")
        settings = entry.copy()
        try:
            length = settings.pop("length")
            rate = settings.pop("rate")
            target = settings.pop("gap_target", None)
            for kind in ("insertion_length", "deletion_length"):
                if kind in settings:
                    distribution = settings[kind]
                    settings[kind] = LengthModel(distribution["name"], tuple(distribution["parameters"]))
            region = Region(length, rate, replace(options, **settings), target)
        except (AttributeError, KeyError, TypeError) as exc:
            raise ValueError(f"invalid region configuration: {entry}") from exc
        regions.append(region)
    return tuple(regions)


def _arguments(parser: argparse.ArgumentParser, argv: list[str] | None) -> argparse.Namespace:
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config", type=Path)
    config, _ = config_parser.parse_known_args(argv)
    if config.config is None:
        return parser.parse_args(argv)
    if argv is None:
        argv = sys.argv[1:]
    try:
        path = config.config
        settings = json.loads(path.read_text(encoding="utf-8-sig"))
        if not isinstance(settings, dict):
            raise ValueError("configuration must be a JSON object")
        allowed = set(vars(parser.parse_args([]))) - {"check", "config"}
        paths = {"regions", "output_dir", "output", "truth", "work_dir", "indelible_bin"}
        lists = {"insertion_length", "deletion_length", "gap_targets"}
        tokens = []
        for key, value in settings.items():
            if key not in allowed:
                raise ValueError(f"unknown configuration option {key!r}")
            if key in paths:
                if not isinstance(value, str) or not value:
                    raise ValueError(f"{key} must be a nonempty path")
                value = str(path.resolve().parent / value)
            if key in lists:
                if not isinstance(value, list) or not value:
                    raise ValueError(f"{key} must be a nonempty array")
                tokens.append("--" + key.replace("_", "-"))
                tokens.extend(map(str, value))
            else:
                if type(value) not in (str, int, float):
                    raise ValueError(f"{key} must be a string or number")
                tokens.append(f"--{key.replace('_', '-')}={value}")
        defaults = parser.parse_args(tokens)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return parser.parse_args(argv, namespace=defaults)


def _check(executable: Path) -> None:
    with TemporaryDirectory(prefix="indelible-check-") as directory:
        result = run_partition(executable, Path(directory) / "run", taxa=("seq1", "seq2"),
                               tree="(seq1:0.1,seq2:0.1);", length=12, seed=13, options=Options())
    print(f"Python {sys.version.split()[0]}, NumPy {np.__version__}")
    print(f"INDELible {result.version}: {executable}")
    print(f"Smoke test passed ({len(result.aligned)} sequences, {len(result.aligned[0])} aligned columns).")


def _single_family(args: argparse.Namespace, executable: Path, options: Options) -> None:
    if args.output_dir is not None:
        if any(path is not None for path in (args.output, args.truth, args.work_dir)):
            raise ValueError("--output-dir cannot be combined with --output, --truth or --work-dir")
        work_dir = args.output_dir
        fasta = work_dir / "alignment.fasta"
        truth = work_dir / "ground_truth.tsv"
    else:
        fasta, truth = args.output, args.truth
        if fasta is None or truth is None:
            raise ValueError("generation requires --output-dir or both --output and --truth")
        work_dir = args.work_dir or fasta.with_suffix(".run")
    raw = work_dir.resolve()
    fasta_path, truth_path = fasta.resolve(), truth.resolve()
    if fasta_path == truth_path or fasta_path in truth_path.parents or truth_path in fasta_path.parents:
        raise ValueError("FASTA and truth must use distinct files")
    for path in (fasta_path, truth_path):
        if path.exists() or path == raw or path in raw.parents:
            raise ValueError(f"expected a new output file: {path}")
        if path.is_relative_to(raw):
            name = path.relative_to(raw).parts[0]
            if name in {"metadata.json", "unaligned.fasta", "calibration"} or name.startswith("region_"):
                raise ValueError(f"output conflicts with a simulator path: {path}")
    if args.replicates is not None:
        raise ValueError("--replicates requires --suite")
    regions = _read_regions(args.regions, options) if args.regions else None
    if regions and args.scenario:
        raise ValueError("choose --regions or --scenario")
    keys = ("size", "length", "divergence", "scenario", "variable_fraction", "motif_length", "seed")
    settings = {key: getattr(args, key) for key in keys if getattr(args, key) is not None}
    if regions:
        settings["scenario"] = "custom"
        settings.setdefault("length", sum(region.length for region in regions))
    simulation = simulate(
        **settings, regions=regions, gaps=args.gaps, gap_by_rate=args.gap_targets, redundancy=args.redundancy,
        executable=executable, work_dir=work_dir, options=options, boundary_width=args.boundary_width,
    )
    fasta.parent.mkdir(parents=True, exist_ok=True)
    truth.parent.mkdir(parents=True, exist_ok=True)
    write_simulation(simulation, fasta, truth, work_dir)
    if reason := simulation.metadata["binary_metrics_reason"]:
        print(reason)
    print(f"Wrote {simulation.alignment.size} sequences and {simulation.alignment.length} aligned columns to {fasta}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    _add_engine_arguments(parser)
    parser.add_argument("--config", type=Path, help="JSON settings with CLI overrides")
    parser.add_argument("--check", action="store_true", help="run a tiny simulator smoke test")
    parser.add_argument("--output-dir", type=Path, help="new directory for FASTA, truth and run metadata")
    parser.add_argument("--output", type=Path, help="write a single family's FASTA to this path")
    parser.add_argument("--truth", type=Path, help="write a single family's truth to this path")
    parser.add_argument("--work-dir", type=Path, help="new raw run directory with --output and --truth")
    design = parser.add_mutually_exclusive_group()
    design.add_argument("--scenario", choices=SCENARIOS + SUPPLEMENTAL_SCENARIOS, help="default: blocks")
    design.add_argument("--regions", type=Path, help="JSON configuration for custom ancestral regions")
    design.add_argument("--suite", choices=SUITES, help="supplemental benchmark family set")
    parser.add_argument("--size", type=int, default=DEFAULT_SIZE, help=f"number of sequences (default: {DEFAULT_SIZE})")
    parser.add_argument("--length", type=int, help=f"ancestral sequence length (default: {DEFAULT_LENGTH})")
    parser.add_argument("--divergence", type=float, help=f"root-to-tip distance (default: {DEFAULT_DIVERGENCE})")
    parser.add_argument("--variable-fraction", type=float)
    parser.add_argument("--motif-length", type=int)
    parser.add_argument("--gaps", type=float, help="target gap occupancy for pilot calibration")
    parser.add_argument("--gap-targets", type=float, nargs=3, metavar=("LOW", "MEDIUM", "HIGH"),
                        help="rate-specific occupancies for three-rate-blocks or d-gaps")
    parser.add_argument("--redundancy", type=float)
    parser.add_argument("--seed", type=int, help="family seed or first suite seed")
    parser.add_argument("--replicates", type=int, help="replicates per suite setting (default: 4)")
    args = _arguments(parser, argv)
    try:
        executable = find_executable(args.indelible_bin)
        if args.check:
            _check(executable)
        elif args.suite:
            if args.output_dir is None or any(value is not None for value in (
                args.output, args.truth, args.work_dir, args.scenario, args.regions, args.divergence,
                args.variable_fraction, args.motif_length, args.gaps, args.redundancy,
            )):
                raise ValueError("suite generation requires --output-dir and uses the suite's region design")
            if args.gap_targets is not None and args.suite != "d-gaps":
                raise ValueError("--gap-targets requires the d-gaps suite or the three-rate-blocks scenario")
            count = generate_suite(
                args.output_dir, args.suite, size=args.size,
                length=DEFAULT_LENGTH if args.length is None else args.length,
                replicates=4 if args.replicates is None else args.replicates, start_seed=args.seed,
                gap_targets=tuple(args.gap_targets) if args.gap_targets is not None else D_GAP_TARGETS,
                executable=executable, options=_engine_options(args), boundary_width=args.boundary_width,
            )
            print(f"Wrote {count} families to {args.output_dir}")
        else:
            _single_family(args, executable, _engine_options(args))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 0



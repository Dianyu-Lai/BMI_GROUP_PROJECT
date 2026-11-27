"""Generate protein MSAs with INDELible and annotate the true alignment."""

import argparse
import csv
from dataclasses import asdict, dataclass, replace
import hashlib
from itertools import groupby
import json
from math import isfinite
from pathlib import Path
from random import Random
from tempfile import TemporaryDirectory

import numpy as np

from msa_conservation.msa import Alignment
from .indelible import LengthModel, Options, find_executable, partition_seed, run_partition


SCENARIOS = ("conserved", "variable", "blocks", "motifs", "gradual",
             "divergent", "indels", "redundant")
SUPPLEMENTAL_SCENARIOS = ("conserved-variable-blocks", "single-motif", "three-rate-blocks")
D_RATES = (0.06, 0.32, 1.0)


@dataclass(frozen=True)
class Simulation:
    alignment: Alignment
    rates: np.ndarray
    truth: np.ndarray
    conserved: np.ndarray
    boundaries: tuple[int, ...]
    metadata: dict | None = None


@dataclass(frozen=True)
class Region:
    length: int
    rate: float
    options: Options | None = None
    gap_target: float | None = None

    def __post_init__(self):
        if type(self.length) is not int or self.length < 1 or not isfinite(self.rate) or self.rate <= 0:
            raise ValueError("each region requires a positive integer length and a finite positive rate")
        if self.gap_target is not None and not 0 <= self.gap_target < 1:
            raise ValueError("region gap targets must be in [0, 1)")


def scenario_rates(scenario: str, length: int, variable_fraction: float, motif_length: int) -> tuple[float, ...]:
    if scenario not in SCENARIOS + SUPPLEMENTAL_SCENARIOS:
        raise ValueError(f"unknown scenario: {scenario}")
    if length < 12 or not 0 < variable_fraction < 1:
        raise ValueError("require at least 12 ancestral positions and variable fraction in (0, 1)")
    low, medium, high = D_RATES
    rates = np.full(length, low)
    if scenario == "variable":
        rates.fill(high)
    elif scenario in ("motifs", "single-motif"):
        if not 1 <= motif_length <= length // 4:
            raise ValueError("motif length must be between 1 and one-quarter of the ancestral length")
        rates.fill(high)
        starts = ([(length - motif_length) // 2] if scenario == "single-motif"
                  else [int(fraction * length) for fraction in (.2, .5, .8)])
        for start in starts:
            rates[start:start + motif_length] = low
    elif scenario == "gradual":
        for bin_indices in np.array_split(np.arange(length), 12):
            phase = 2 * np.pi * float(bin_indices.mean()) / (length - 1)
            rates[bin_indices] = round(low + (high - low) * (1 - np.cos(phase)) / 2, 12)
    elif scenario != "conserved":
        block = min(length - 6, max(2, round(length * variable_fraction)))
        left = (length - block) // 2
        right = left + block
        if scenario == "three-rate-blocks":
            medium_width = left // 2
            rates[left - medium_width:left] = medium
            rates[right:right + medium_width] = medium
        elif scenario != "conserved-variable-blocks":
            rates[max(2, left // 2):left] = medium
        rates[left:right] = high
    return tuple(rates.tolist())


def _tree(taxa: tuple[str, ...], seed: int, redundancy: float) -> tuple:
    """Build a rooted, ultrametric binary topology with optional near-tip relatives."""
    order = list(taxa)
    Random(seed).shuffle(order)
    cluster_size = min(len(order) - 1, round(redundancy * len(order)) + 1)
    cluster = frozenset(order[:cluster_size]) if redundancy and cluster_size >= 2 else frozenset()

    def split(names: list[str], parent: float) -> tuple:
        if len(names) == 1:
            return (names[0], 1 - parent, ())
        depth = parent + (1 - parent) * (.98 if frozenset(names) == cluster else .4)
        middle = len(names) // 2
        return ("", depth - parent, (split(names[:middle], depth), split(names[middle:], depth)))

    middle = len(cluster) or len(order) // 2
    return ("", 0., (split(order[:middle], 0.), split(order[middle:], 0.)))


def _newick(node: tuple, scale: float) -> str:
    name, length, children = node
    if children:
        name = "(" + ",".join(_newick(child, scale) for child in children) + ")"
    return f"{name}:{length * scale:.15g}" if length else name + ";"


def _gap_statistics(gaps: np.ndarray) -> dict:
    gap_runs: dict[int, int] = {}
    for row in gaps:
        for gap, values in groupby(row):
            if gap:
                length = sum(1 for _ in values)
                gap_runs[length] = gap_runs.get(length, 0) + 1
    return {"gap_occupancy": float(gaps.mean()),
            "per_sequence_gap_occupancy": gaps.mean(axis=1).tolist(), "gap_run_lengths": gap_runs}


def _matrix(aligned: tuple[str, ...]) -> np.ndarray:
    return np.frombuffer("".join(aligned).encode("ascii"), dtype="S1").reshape(len(aligned), -1)


def _statistics(aligned: tuple[str, ...]) -> dict:
    matrix = _matrix(aligned)
    gaps = matrix == b"-"
    identities = []
    for index, row in enumerate(matrix[:-1]):
        compared = ~(gaps[index] | gaps[index + 1:])
        totals = compared.sum(axis=1)
        matches = ((row == matrix[index + 1:]) & compared).sum(axis=1)
        valid = totals > 0
        identities.extend((matches[valid] / totals[valid]).tolist())
    identity = float(np.mean(identities)) if identities else None
    return {"mean_pairwise_identity": identity,
            "mean_pairwise_p_distance": None if identity is None else 1 - identity,
            "comparable_pairs": len(identities), **_gap_statistics(gaps)}


def _calibrate(executable: Path, directory: Path, *, taxa: tuple[str, ...], tree: str,
               length: int, seed: int, region: int, target: float, options: Options) -> tuple[Options, list[dict]]:
    if not 0 <= target < 1:
        raise ValueError("target gap occupancy must be in [0, 1)")
    if target == 0:
        return replace(options, insertion_rate=0., deletion_rate=0.), []
    total_rate = options.insertion_rate + options.deletion_rate
    insert_share = options.insertion_rate / total_rate if total_rate else .5
    delete_share = 1 - insert_share
    lower, upper = 0., 1.5
    best = replace(options, insertion_rate=0., deletion_rate=0.)
    best_error = target
    trials = []
    for trial in range(7):
        candidate = (lower + upper) / 2
        settings = replace(options, insertion_rate=candidate * insert_share,
                           deletion_rate=candidate * delete_share)
        observations = []
        for pilot in range(1, 6):
            run = run_partition(executable, directory / f"trial_{trial + 1}_{pilot}", taxa=taxa,
                                tree=tree, length=length, seed=partition_seed(seed, region, pilot), options=settings)
            stats = _gap_statistics(_matrix(run.aligned) == b"-")
            observations.append({"gap_occupancy": stats["gap_occupancy"],
                                 "gap_run_lengths": stats["gap_run_lengths"], "indel_events": run.indel_events,
                                 "seconds": run.seconds, "exit_status": 0})
        measured = float(np.mean([item["gap_occupancy"] for item in observations]))
        trials.append({"insertion_rate": settings.insertion_rate, "deletion_rate": settings.deletion_rate,
                       "mean_gap_occupancy": measured, "pilots": observations})
        if abs(measured - target) < best_error:
            best, best_error = settings, abs(measured - target)
        if measured < target:
            lower = candidate
        else:
            upper = candidate
        if best_error <= .01:
            break
    return best, trials


def simulate(*, size: int = 50, length: int = 240, divergence: float = .7,
             scenario: str = "blocks", variable_fraction: float = .4, motif_length: int = 5,
             gaps: float | None = None, redundancy: float | None = None,
             gap_by_rate: tuple[float, float, float] | None = None, seed: int = 1,
             executable: Path | None = None, work_dir: Path | None = None,
             options: Options = Options(), boundary_width: int = 0,
             ancestral_rates: tuple[float, ...] | None = None,
             regions: tuple[Region, ...] | None = None) -> Simulation:
    """Simulate each designed region on the same topology with scaled branches."""
    if size < 2 or not isfinite(divergence) or divergence < 0 or boundary_width < 0:
        raise ValueError("require at least two sequences and nonnegative divergence and boundary width")
    redundancy = (.7 if scenario == "redundant" else 0.) if redundancy is None else redundancy
    if not 0 <= redundancy < 1:
        raise ValueError("redundancy must be in [0, 1)")
    if gaps is None and scenario == "indels" and not (options.insertion_rate or options.deletion_rate):
        gaps = .2
    if gap_by_rate is not None and (gaps is not None or scenario != "three-rate-blocks"
                                   or len(gap_by_rate) != len(D_RATES)):
        raise ValueError("rate-specific gap targets require three-rate-blocks without a global gap target")
    if regions is None:
        rates = (scenario_rates(scenario, length, variable_fraction, motif_length)
                 if ancestral_rates is None else ancestral_rates)
        regions = tuple(Region(sum(1 for _ in positions), rate) for rate, positions in groupby(rates))
    if not regions or sum(region.length for region in regions) != length:
        raise ValueError("region lengths must sum to the ancestral sequence length")
    executable = find_executable(executable)
    if work_dir is None:
        with TemporaryDirectory(prefix="bmi3-indelible-") as directory:
            return _simulate(executable, Path(directory) / "raw", size, divergence, scenario, redundancy,
                             gaps, gap_by_rate, seed, options, boundary_width, regions)
    return _simulate(executable, work_dir, size, divergence, scenario, redundancy,
                     gaps, gap_by_rate, seed, options, boundary_width, regions)


def _simulate(executable, work_dir, size, divergence, scenario, redundancy, gaps,
              gap_by_rate, seed, options, boundary_width, regions) -> Simulation:
    work_dir.mkdir(parents=True)
    taxa = tuple(f"seq{i + 1}" for i in range(size))
    topology = _tree(taxa, seed, redundancy)
    divergence *= 2.5 if scenario == "divergent" else 1.
    targets = dict(zip(D_RATES, gap_by_rate)) if gap_by_rate is not None else {}
    segments, partitions, final_rates = [], [], []
    calibrations = {}
    root_start = 1
    for region, design in enumerate(regions, 1):
        root_length, rate = design.length, design.rate
        tree = _newick(topology, divergence * rate)
        target = targets.get(rate, gaps) if design.gap_target is None else design.gap_target
        settings, trials = design.options or options, []
        name = f"region_{region:02d}"
        calibration_path = None
        if target is not None:
            key = (root_length, rate, target, settings)
            if key not in calibrations:
                calibration_path = Path("calibration") / name
                calibrated, trials = _calibrate(executable, work_dir / calibration_path, taxa=taxa,
                                                tree=tree, length=root_length, seed=seed, region=region,
                                                target=target, options=settings)
                calibrations[key] = calibrated, trials, calibration_path
            settings, trials, calibration_path = calibrations[key]
        engine_seed = partition_seed(seed, region)
        segment = run_partition(executable, work_dir / name, taxa=taxa, tree=tree,
                                length=root_length, seed=engine_seed, options=settings)
        rows = tuple(row.upper() for row in segment.aligned)
        aligned_length = len(rows[0])
        aligned_start = len(final_rates) + 1
        final_rates.extend([rate] * aligned_length)
        segments.append(rows)
        partitions.append({"id": region, "ancestral_start": root_start, "ancestral_length": root_length,
                           "aligned_start": aligned_start, "aligned_length": aligned_length, "rate": rate,
                           "tree": tree, "engine_seed": engine_seed, "options": asdict(settings),
                           "gap_target": target, "calibration": trials, "statistics": _statistics(rows),
                           "calibration_path": calibration_path.as_posix() if calibration_path and trials else None,
                           "indel_events": segment.indel_events,
                           "control_file": f"{name}/control.txt", "exit_status": 0, "seconds": segment.seconds})
        root_start += root_length
    aligned = tuple("".join(parts) for parts in zip(*segments))
    if any(not row.strip("-") for row in aligned):
        raise ValueError(f"a combined descendant has no residues in {work_dir}")
    rates = np.array(final_rates)
    conserved = rates <= .18
    boundaries = tuple((np.flatnonzero(np.diff(rates)) + 1).tolist())
    metadata = {"simulator": "INDELible", "version": segment.version,
                "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
                "seed": seed, "scenario": scenario,
                "sequence_count": size, "ancestral_length": root_start - 1, "aligned_length": len(rates),
                "divergence": divergence, "tree": _newick(topology, divergence), "taxa": taxa,
                "redundancy": redundancy, "boundary_width": boundary_width, "partitions": partitions,
                "statistics": _statistics(aligned), "simulation_seconds": sum(p["seconds"] for p in partitions),
                "binary_metrics_reason": (None if conserved.any() and not conserved.all()
                                          else "AUROC/AUPRC undefined because designed truth contains one class")}
    return Simulation(Alignment(taxa, aligned), rates, np.exp(-divergence * rates), conserved, boundaries, metadata)


def write_simulation(simulation: Simulation, fasta: Path, truth: Path, work_dir: Path) -> None:
    metadata = simulation.metadata
    alignment = simulation.alignment
    boundaries = set(simulation.boundaries)
    width = metadata["boundary_width"]
    joins = [part["aligned_start"] - 1 for part in metadata["partitions"][1:]]
    mask = {index for join in joins for index in range(join - width, join + width)}
    with fasta.open("x", encoding="ascii", newline="\n") as handle:
        for name, row in zip(alignment.ids, alignment.sequences):
            handle.write(f">{name}\n{row}\n")
    with (work_dir / "unaligned.fasta").open("x", encoding="ascii", newline="\n") as handle:
        for name, row in zip(alignment.ids, alignment.sequences):
            handle.write(f">{name}\n{row.replace('-', '')}\n")
    with truth.open("x", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("position", "rate", "conservation", "conserved", "boundary_before", "partition_id",
                         "region_id", "mutation_rate_class", "boundary_mask", "observed_gap_fraction"))
        for part in metadata["partitions"]:
            start = part["aligned_start"] - 1
            for index in range(start, start + part["aligned_length"]):
                rate_class = "low" if simulation.conserved[index] else "high" if simulation.rates[index] >= 1 else "medium"
                gap_fraction = sum(row[index] == "-" for row in alignment.sequences) / alignment.size
                writer.writerow((index + 1, simulation.rates[index], simulation.truth[index],
                                 int(simulation.conserved[index]), int(index in boundaries), part["id"],
                                 part["id"], rate_class, int(index in mask), gap_fraction))
    (work_dir / "metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def add_engine_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--indelible-bin", type=Path)
    parser.add_argument("--model", choices=("WAG", "LG", "JTT"), default="WAG")
    parser.add_argument("--method", type=int, choices=(1, 2), default=1)
    parser.add_argument("--insertion-rate", type=float, default=0.)
    parser.add_argument("--deletion-rate", type=float, default=0.)
    parser.add_argument("--insertion-length", nargs="+", default=("POW", "1.7", "15"), metavar="MODEL_OR_PARAMETER")
    parser.add_argument("--deletion-length", nargs="+", default=("POW", "1.7", "15"), metavar="MODEL_OR_PARAMETER")
    parser.add_argument("--timeout", type=float, default=60., help="seconds per simulator process")
    parser.add_argument("--boundary-width", type=int, default=0, help="mask this many columns on each side")


def engine_options(args: argparse.Namespace) -> Options:
    return Options(args.model, args.method, args.insertion_rate, args.deletion_rate,
                   LengthModel(args.insertion_length[0], tuple(map(float, args.insertion_length[1:]))),
                   LengthModel(args.deletion_length[0], tuple(map(float, args.deletion_length[1:]))), args.timeout)


def read_regions(path: Path, options: Options) -> tuple[Region, ...]:
    entries = json.loads(path.read_text(encoding="utf-8"))
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
        except (KeyError, TypeError) as exc:
            raise ValueError(f"invalid region configuration: {entry}") from exc
        regions.append(region)
    return tuple(regions)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_engine_arguments(parser)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, help="new path for raw runs and metadata")
    parser.add_argument("--scenario", choices=SCENARIOS + SUPPLEMENTAL_SCENARIOS, default="blocks")
    parser.add_argument("--regions", type=Path, help="JSON configuration for custom ancestral regions")
    parser.add_argument("--size", type=int, default=50)
    parser.add_argument("--length", type=int, help="ancestral sequence length (default: 240)")
    parser.add_argument("--divergence", type=float, default=.7)
    parser.add_argument("--variable-fraction", type=float, default=.4)
    parser.add_argument("--motif-length", type=int, default=5)
    parser.add_argument("--gaps", type=float, help="target gap occupancy for pilot calibration")
    parser.add_argument("--redundancy", type=float)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    work_dir = args.work_dir or args.output.with_suffix(".run")
    try:
        if args.output.resolve() == args.truth.resolve() or args.output.exists() or args.truth.exists():
            raise ValueError("FASTA and truth must use distinct new files")
        options = engine_options(args)
        regions = read_regions(args.regions, options) if args.regions else None
        length = args.length if args.length is not None else (sum(region.length for region in regions) if regions else 240)
        simulator = simulate(scenario="custom" if regions else args.scenario, size=args.size,
                             length=length, divergence=args.divergence, regions=regions,
                             variable_fraction=args.variable_fraction, motif_length=args.motif_length, gaps=args.gaps,
                             redundancy=args.redundancy, seed=args.seed, executable=args.indelible_bin, work_dir=work_dir,
                             options=options, boundary_width=args.boundary_width)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.truth.parent.mkdir(parents=True, exist_ok=True)
        write_simulation(simulator, args.output, args.truth, work_dir)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    if reason := simulator.metadata["binary_metrics_reason"]:
        print(reason)
    print(f"Wrote {simulator.alignment.size} sequences and {simulator.alignment.length} aligned columns to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

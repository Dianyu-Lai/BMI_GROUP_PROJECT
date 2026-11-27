"""Generate independent synthetic protein families with known mutation-rate profiles."""

from dataclasses import dataclass
import argparse
import csv
from pathlib import Path

import numpy as np

from msa_conservation.msa import AMINO_ACIDS, Alignment


SCENARIOS = ("conserved", "variable", "blocks", "motifs", "gradual",
             "divergent", "indels", "redundant")
SUPPLEMENTAL_SCENARIOS = ("conserved-variable-blocks", "single-motif", "three-rate-blocks")
CHEMICAL_CLASSES = ("AVLIM", "FWY", "STNQ", "KR", "DE", "C", "G", "H", "P")
MUTATIONS = ("mixed", "conservative", "radical")
D_RATES = (0.06, 0.32, 1.0)


@dataclass(frozen=True)
class Simulation:
    alignment: Alignment
    rates: np.ndarray
    truth: np.ndarray
    conserved: np.ndarray
    boundaries: tuple[int, ...]
    gap_probabilities: np.ndarray | None = None


def simulate(
    *, size: int = 50, length: int = 240, divergence: float = 0.7,
    scenario: str = "blocks", variable_fraction: float = 0.4,
    gaps: float | None = None, redundancy: float | None = None,
    gap_by_rate: tuple[float, float, float] | None = None, gap_seed: int | None = None,
    mutation: str = "mixed", motif_length: int = 5, seed: int = 1,
) -> Simulation:
    """Sample from one ancestor with independent mutation masks and optional short indels."""
    if size < 2 or length < 12 or not np.isfinite(divergence) or divergence <= 0:
        raise ValueError("require at least two sequences, 12 sites and positive finite divergence")
    if scenario not in SCENARIOS + SUPPLEMENTAL_SCENARIOS or mutation not in MUTATIONS:
        raise ValueError("unknown scenario or substitution regime")
    if not 0 < variable_fraction < 1 or (
        scenario in ("motifs", "single-motif") and not 1 <= motif_length <= length // 4
    ):
        raise ValueError("invalid variable fraction or motif length")
    gaps = (0.2 if scenario == "indels" else 0.0) if gaps is None else gaps
    redundancy = (0.7 if scenario == "redundant" else 0.0) if redundancy is None else redundancy
    if not (0 <= gaps < 1 and 0 <= redundancy < 1):
        raise ValueError("gap and redundancy fractions must be in [0, 1)")
    if gap_by_rate is not None and (
        scenario != "three-rate-blocks" or gaps or len(gap_by_rate) != len(D_RATES)
        or not all(0 <= probability < 1 for probability in gap_by_rate)
    ):
        raise ValueError("rate-specific gaps require three probabilities in [0, 1) and the three-rate scenario without global gaps")
    if gap_seed is not None and gap_by_rate is None:
        raise ValueError("gap_seed requires gap_by_rate")
    effective_divergence = divergence * (2.5 if scenario == "divergent" else 1)
    boundaries = ()
    if scenario in ("conserved", "variable"):
        rates = np.full(length, .06 if scenario == "conserved" else 1.)
    elif scenario in ("motifs", "single-motif"):
        rates = np.ones(length)
        if scenario == "single-motif":
            starts = [(length - motif_length) // 2]
        else:
            starts = [int(fraction * length) for fraction in (.2, .5, .8)]
        for start in starts:
            rates[start:min(start + motif_length, length)] = .06
        boundaries = tuple(np.flatnonzero(np.diff(rates)) + 1)
    elif scenario == "gradual":
        rates = .06 + .94 * (1 - np.cos(np.linspace(0, 2 * np.pi, length))) / 2
    else:
        rates = np.full(length, .06)
        block = min(length - 6, max(2, round(length * variable_fraction)))
        left = (length - block) // 2
        if scenario == "three-rate-blocks":
            medium_width = left // 2
            right = left + block
            rates[left - medium_width:left] = .32
            rates[right:right + medium_width] = .32
        elif scenario != "conserved-variable-blocks":
            moderate = max(2, left // 2)
            rates[moderate:left] = .32
        rates[left:left + block] = 1.
        boundaries = tuple(np.flatnonzero(np.diff(rates)) + 1)

    rng = np.random.default_rng(seed)
    alphabet = np.array(list(AMINO_ACIDS))
    ancestor = rng.integers(len(alphabet), size=length)
    unique = max(2, round(size * (1 - redundancy)))
    scores = effective_divergence * rates
    truth = np.exp(-scores)
    changed = rng.random((unique, length)) < -np.expm1(-scores)
    descendants = np.broadcast_to(ancestor, changed.shape).copy()
    for index, residue in enumerate(AMINO_ACIDS):
        cluster = next(group for group in CHEMICAL_CLASSES if residue in group)
        alternatives = [j for j, aa in enumerate(AMINO_ACIDS) if aa != residue and (
            mutation == "mixed" or (mutation == "conservative" and aa in cluster)
            or (mutation == "radical" and aa not in cluster)
        )]
        # Singleton classes keep the site mutation probability through uniform fallback.
        if not alternatives:
            alternatives = [j for j in range(len(alphabet)) if j != index]
        mask = changed & (ancestor == index)
        descendants[mask] = rng.choice(alternatives, size=int(mask.sum()))

    # A shared subset of positions models aligned inserts absent from many rows.
    # Per-row short deletion runs supply the remaining target coverage loss.
    chars = alphabet[descendants]
    if gaps:
        target = round(gaps * length)
        insert_sites = rng.choice(length, size=target // 3, replace=False)
        chars[:, insert_sites] = np.where(
            rng.random((unique, len(insert_sites))) < .75, "-", chars[:, insert_sites]
        )
        for row in chars:
            missing = int(np.count_nonzero(row == "-"))
            while missing < target:
                start = int(rng.integers(length))
                stop = min(length, start + min(target - missing, int(rng.integers(2, 7))))
                missing += int(np.count_nonzero(row[start:stop] != "-"))
                row[start:stop] = "-"
    gap_probabilities = None
    if gap_by_rate is not None:
        probability_by_rate = dict(zip(D_RATES, gap_by_rate))
        gap_probabilities = np.array([probability_by_rate[rate] for rate in rates])
        gap_rng = rng if gap_seed is None else np.random.default_rng(gap_seed)
        chars[gap_rng.random(chars.shape) < gap_probabilities] = "-"
    if unique < size:
        chars = np.concatenate((chars, chars[rng.integers(unique, size=size - unique)]))
    alignment = Alignment(tuple(f"seq{i+1}" for i in range(size)),
                          tuple("".join(row) for row in chars))
    return Simulation(alignment, rates, truth, rates <= .18, boundaries, gap_probabilities)


def write_simulation(simulation: Simulation, fasta_path: Path, truth_path: Path) -> None:
    """Export the aligned FASTA and one row of ground truth per column."""
    with fasta_path.open("w", encoding="utf-8") as handle:
        for name, seq in zip(simulation.alignment.ids, simulation.alignment.sequences):
            handle.write(f">{name}\n{seq}\n")
    with truth_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        fields = ("position", "rate", "conservation", "conserved", "boundary_before")
        gap_probabilities = simulation.gap_probabilities
        if gap_probabilities is not None:
            fields += ("gap_probability", "observed_gap_fraction")
            gap_fractions = [column.count("-") / simulation.alignment.size
                             for column in zip(*simulation.alignment.sequences)]
        writer.writerow(fields)
        boundaries = set(simulation.boundaries)
        for i, (rate, target, conserved) in enumerate(zip(
            simulation.rates, simulation.truth, simulation.conserved
        )):
            row = (i + 1, rate, target, int(conserved), int(i in boundaries))
            if gap_probabilities is not None:
                row += (gap_probabilities[i], gap_fractions[i])
            writer.writerow(row)


def main() -> int:
    """Write a simulated aligned FASTA file and its exact truth on aligned positions."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--scenario", choices=SCENARIOS + SUPPLEMENTAL_SCENARIOS, default="blocks")
    parser.add_argument("--mutation", choices=MUTATIONS, default="mixed")
    parser.add_argument("--size", type=int, default=50)
    parser.add_argument("--length", type=int, default=240)
    parser.add_argument("--divergence", type=float, default=.7)
    parser.add_argument("--variable-fraction", type=float, default=.4)
    parser.add_argument("--motif-length", type=int, default=5)
    parser.add_argument("--gaps", type=float)
    parser.add_argument("--redundancy", type=float)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    options = vars(args).copy()
    output, truth = options.pop("output"), options.pop("truth")
    if output.resolve() == truth.resolve():
        parser.error("FASTA and truth output paths must be distinct")
    try:
        write_simulation(simulate(**options), output, truth)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

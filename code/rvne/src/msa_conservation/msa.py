"""Aligned FASTA parsing, sequence weights, and weighted column frequencies."""

from collections import Counter
from dataclasses import dataclass
from math import fsum, isfinite
from typing import Iterable, Iterator, Sequence


AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"
AA_INDEX = {residue: index for index, residue in enumerate(AMINO_ACIDS)}
GAPS = frozenset("-.")
ALLOWED = frozenset(AMINO_ACIDS + "BJOUXZ-.")
ENCODING = {**AA_INDEX, **dict.fromkeys("BJOUXZ", -1), **dict.fromkeys(GAPS, -2)}


@dataclass(frozen=True)
class Alignment:
    ids: tuple[str, ...]
    sequences: tuple[str, ...]

    def __post_init__(self) -> None:
        ids = tuple(self.ids)
        sequences = tuple(sequence.upper() for sequence in self.sequences)
        if not sequences or not sequences[0]:
            raise ValueError("alignment must contain at least one nonempty sequence")
        if len(ids) != len(sequences) or any(not name.strip() for name in ids):
            raise ValueError("every sequence must have a nonempty FASTA identifier")
        if len(set(ids)) != len(ids):
            raise ValueError("FASTA identifiers must be unique")

        length = len(sequences[0])
        for name, sequence in zip(ids, sequences):
            if len(sequence) != length:
                raise ValueError(f"{name!r} has length {len(sequence)}, expected {length}")
            invalid = sorted(set(sequence) - ALLOWED)
            if invalid:
                raise ValueError(f"{name!r} contains unsupported symbols: {', '.join(invalid)}")
        object.__setattr__(self, "ids", ids)
        object.__setattr__(self, "sequences", sequences)

    @property
    def length(self) -> int:
        return len(self.sequences[0])

    @property
    def size(self) -> int:
        return len(self.sequences)


@dataclass(frozen=True)
class ColumnProfile:
    frequencies: tuple[float, ...]
    nongap_coverage: float
    canonical_coverage: float


def read_fasta_records(lines: Iterable[str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Read FASTA records without requiring equal sequence lengths."""
    ids: list[str] = []
    sequences: list[str] = []
    fragments: list[str] = []

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if ids:
                sequences.append("".join(fragments))
            description = line[1:].split()
            if not description:
                raise ValueError(f"line {line_number}: empty FASTA header")
            ids.append(description[0])
            fragments = []
        else:
            if not ids:
                raise ValueError(f"line {line_number}: sequence before FASTA header")
            fragments.append("".join(line.split()))

    if ids:
        sequences.append("".join(fragments))
    return tuple(ids), tuple(sequences)


def read_fasta(lines: Iterable[str]) -> Alignment:
    return Alignment(*read_fasta_records(lines))


def henikoff_weights(alignment: Alignment) -> tuple[float, ...]:
    """Sum 1 / (distinct residues × copies of this residue) per valid site."""
    weights = [0.0] * alignment.size
    for column in zip(*alignment.sequences):
        counts = Counter(residue for residue in column if residue in AA_INDEX)
        diversity = len(counts)
        for i, residue in enumerate(column):
            if residue in counts:
                weights[i] += 1.0 / (diversity * counts[residue])
    total = fsum(weights)
    if total == 0:
        return (1.0 / alignment.size,) * alignment.size
    return tuple(weight / total for weight in weights)


def column_profiles(
    alignment: Alignment, weights: Sequence[float] | None = None
) -> Iterator[ColumnProfile]:
    """Yield frequencies conditional on observing a canonical amino acid."""
    selected = tuple(weights) if weights is not None else (1.0,) * alignment.size
    if len(selected) != alignment.size:
        raise ValueError("provide exactly one weight per aligned sequence")
    if any(not isfinite(weight) or weight < 0 for weight in selected):
        raise ValueError("sequence weights must be finite and nonnegative")
    scale = max(selected)
    if scale == 0:
        raise ValueError("at least one sequence weight must be positive")
    normalized = tuple(weight / scale for weight in selected)
    total_weight = fsum(normalized)

    for column in zip(*alignment.sequences):
        counts = [0.0] * len(AMINO_ACIDS)
        nongap_weights = []
        for residue, weight in zip(column, normalized):
            index = ENCODING[residue]
            if index == -2:
                continue
            nongap_weights.append(weight)
            if index >= 0:
                counts[index] += weight

        canonical_weight = fsum(counts)
        frequencies = tuple(value / canonical_weight for value in counts) if canonical_weight else tuple(counts)
        yield ColumnProfile(
            frequencies=frequencies,
            nongap_coverage=min(1.0, fsum(nongap_weights) / total_weight),
            canonical_coverage=min(1.0, canonical_weight / total_weight),
        )

"""Run INDELible 1.03 and validate its raw partition outputs."""

import hashlib
from dataclasses import dataclass
from math import isfinite
import os
from pathlib import Path
import re
import shutil
import subprocess
from time import perf_counter


RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class LengthModel:
    name: str = "POW"
    parameters: tuple[float, ...] = (1.7, 15)

    def __post_init__(self):
        if self.name not in {"NB", "POW", "LAV", "USER"}:
            raise ValueError("indel-length model must be NB, POW, LAV or USER")
        args = self.parameters
        if not args or not all(isfinite(value) for value in args):
            raise ValueError(f"{self.name} requires finite parameters")
        if self.name == "NB":
            valid = len(args) == 2 and 0 <= args[0] < 1 and args[1] >= 1 and args[1] % 1 == 0
        elif self.name == "USER":
            valid = min(args) >= 0 and isfinite(sum(args)) and sum(args) > 0
        else:
            valid = args[0] > 1 and (
                (len(args) == 1 and self.name == "POW")
                or (len(args) == 2 and args[1] > 1 and args[1] % 1 == 0)
            )
        if not valid:
            raise ValueError(f"invalid parameters for {self.name}: {args}")

    def stage(self, directory: Path, kind: str) -> str:
        if self.name == "USER":
            filename = f"{kind}_lengths.txt"
            (directory / filename).write_text(" ".join(map(str, self.parameters)) + "\n", encoding="ascii")
            return f"USER {filename}"
        return " ".join((self.name, *(f"{value:g}" for value in self.parameters)))


@dataclass(frozen=True)
class Options:
    model: str = "WAG"
    method: int = 1
    insertion_rate: float = 0.0
    deletion_rate: float = 0.0
    insertion_length: LengthModel = LengthModel()
    deletion_length: LengthModel = LengthModel()
    timeout: float = 60.0

    def __post_init__(self):
        if self.model not in {"WAG", "LG", "JTT"} or self.method not in {1, 2}:
            raise ValueError("choose WAG, LG or JTT and Method 1 or 2")
        if not all(isfinite(rate) and rate >= 0 for rate in (self.insertion_rate, self.deletion_rate)):
            raise ValueError("insertion and deletion rates must be finite and nonnegative")
        if not isfinite(self.timeout) or self.timeout <= 0:
            raise ValueError("timeout must be finite and positive")


@dataclass(frozen=True)
class PartitionOutput:
    aligned: tuple[str, ...]
    unaligned: tuple[str, ...]
    version: str
    seconds: float
    indel_events: dict


def find_executable(path: str | Path | None = None) -> Path:
    configured = str(path) if path is not None else os.environ.get("INDELIBLE_BIN")
    candidate = Path(configured).expanduser() if configured else None
    executable = str(candidate) if candidate and candidate.is_file() else shutil.which(configured or "indelible")
    if not executable:
        detail = f" at {configured!r}" if configured else " on PATH"
        raise FileNotFoundError(f"INDELible 1.03 executable not found{detail}. Set INDELIBLE_BIN or pass --indelible-bin.")
    return Path(executable).resolve()


def partition_seed(seed: int, partition: int, pilot: int = 0) -> int:
    digest = hashlib.sha256(f"{seed}/{partition}/{pilot}".encode("ascii")).digest()
    return 1 + int.from_bytes(digest[:8], "big") % (2**31 - 2)


def _read_fasta(path: Path, taxa: tuple[str, ...], alphabet: frozenset[str]) -> tuple[str, ...]:
    records: dict[str, list[str]] = {}
    name = None
    for line in path.read_text(encoding="ascii").splitlines():
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].strip().split(maxsplit=1)
            if not name or name[0] in records:
                raise ValueError(f"empty or duplicate FASTA identifier in {path}")
            name = name[0]
            records[name] = []
        elif line:
            if name is None or not set(line.upper()) <= alphabet:
                raise ValueError(f"invalid FASTA sequence in {path}")
            records[name].append(line)
    if set(records) != set(taxa):
        raise ValueError(f"taxon identifiers do not match the guide tree in {path}")
    return tuple("".join(records[taxon]) for taxon in taxa)


def run_partition(executable: Path, directory: Path, *, taxa: tuple[str, ...],
                  tree: str, length: int, seed: int, options: Options) -> PartitionOutput:
    directory.mkdir(parents=True)
    insertion = options.insertion_length.stage(directory, "insertion")
    deletion = options.deletion_length.stage(directory, "deletion")
    control = f"""[TYPE] AMINOACID {options.method}
[SETTINGS]
[output] FASTA
[randomseed] {seed}
[fileperrep] TRUE
[ancestralprint] NEW
[printrates] TRUE
[insertaslowercase] TRUE
[markdeletedinsertions] FALSE
[MODEL] protein
[submodel] {options.model}
[insertmodel] {insertion}
[deletemodel] {deletion}
[insertrate] {options.insertion_rate:.15g}
[deleterate] {options.deletion_rate:.15g}
[TREE] guide {tree}
[PARTITIONS] region [guide protein {length}]
[EVOLVE] region 1 simulation
"""
    (directory / "control.txt").write_text(control, encoding="ascii")
    staged = directory / ("indelible.exe" if os.name == "nt" else "indelible")
    shutil.copy2(executable, staged)
    start = perf_counter()
    try:
        with (directory / "stdout.log").open("wb") as stdout, (directory / "stderr.log").open("wb") as stderr:
            result = subprocess.run([str(staged.resolve())], cwd=directory, stdin=subprocess.DEVNULL,
                                    stdout=stdout, stderr=stderr, timeout=options.timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        raise ValueError(f"INDELible timed out after {options.timeout:g}s. Logs and control file: {directory}") from exc
    except OSError as exc:
        raise ValueError(f"could not execute INDELible at {executable}: {exc}. Run files: {directory}") from exc
    finally:
        staged.unlink()
    elapsed = perf_counter() - start
    stdout = (directory / "stdout.log").read_text(encoding="ascii", errors="replace")
    if result.returncode:
        raise ValueError(f"INDELible exited with status {result.returncode}. Logs and control file: {directory}")
    version = re.search(r"INDELible\s+V([\d.]+)", stdout, re.IGNORECASE)
    if not version or version.group(1) != "1.03":
        raise ValueError(f"expected INDELible 1.03. Check stdout.log in {directory}")
    aligned = _read_fasta(directory / "simulation_TRUE_1.fas", taxa, RESIDUES | {"-"})
    unaligned = _read_fasta(directory / "simulation_1.fas", taxa, RESIDUES)
    if not aligned[0] or len({len(row) for row in aligned}) != 1 or not any(unaligned):
        raise ValueError(f"empty or inconsistent partition alignment in {directory}")
    if any(row.replace("-", "").upper() != raw.upper() for row, raw in zip(aligned, unaligned)):
        raise ValueError(f"true alignment and unaligned output disagree in {directory}")
    log = (directory / "LOG.txt").read_text(encoding="ascii", errors="replace")
    events = {}
    for kind in ("insertion", "deletion"):
        count = re.search(rf"Number of {kind} events\s+(\d+)", log)
        average = re.search(rf"Actual average {kind} length\s+(\S+)", log)
        events[kind] = {"count": int(count.group(1)) if count else 0,
                        "mean_event_length": float(average.group(1))
                        if count and int(count.group(1)) and average else None}
    return PartitionOutput(aligned, unaligned, version.group(1), elapsed, events)

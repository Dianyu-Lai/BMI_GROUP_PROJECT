"""Align an unaligned protein family and compare conservation in reference coordinates."""

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
from statistics import StatisticsError, correlation
import subprocess
from tempfile import TemporaryDirectory
from textwrap import wrap
from time import perf_counter
from typing import Iterable

from . import rvne_cli
from .msa import AMINO_ACIDS, GAPS, Alignment, read_fasta, read_fasta_records


MAFFT_MODES = {
    "auto": ("--auto",),
    "linsi": ("--localpair", "--maxiterate", "1000"),
    "ginsi": ("--globalpair", "--maxiterate", "1000"),
    "einsi": ("--genafpair", "--maxiterate", "1000", "--ep", "0"),
    "fftns2": ("--retree", "2", "--maxiterate", "0"),
}
EXECUTABLES = {"mafft": ("mafft",), "clustalw": ("clustalw2", "clustalw")}
VERSION_PATTERNS = {
    "mafft": re.compile(r"\bv(\d+\.\d+(?:\.\d+)?)\b", re.IGNORECASE),
    "clustalw": re.compile(r"\bCLUSTAL(?:\s*W)?\s+\(?(?:version\s+)?(\d+\.\d+(?:\.\d+)?)\b", re.IGNORECASE),
}


def read_proteins(data: bytes) -> dict[str, str]:
    ids, sequences = read_fasta_records(data.decode("utf-8-sig").splitlines())
    if len(ids) < 2 or len(ids) != len(set(ids)):
        raise ValueError("input FASTA requires at least two sequences with unique IDs")
    proteins = {}
    for name, sequence in zip(ids, sequences):
        if not sequence:
            raise ValueError(f"{name!r} has an empty sequence")
        if set(sequence) & GAPS:
            raise ValueError(f"{name!r} contains gaps. Provide unaligned sequences")
        residues = sequence.upper()
        invalid = sorted(set(residues) - set(AMINO_ACIDS))
        if invalid or not sequence.isascii():
            raise ValueError(f"{name!r} must use the 20 canonical amino acids (found {', '.join(invalid) or 'non-ASCII characters'})")
        proteins[name] = residues
    return proteins


def write_fasta(path: Path, ids: Iterable[str], sequences: Iterable[str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in zip(ids, sequences, strict=True):
            handle.write(f">{name}\n")
            handle.write("\n".join(wrap(sequence, 80)) + "\n")


def resolve_executable(aligner: str, override: str | None) -> Path:
    for name in (override,) if override else EXECUTABLES[aligner]:
        executable = shutil.which(name)
        if executable:
            return Path(executable).absolute()
    raise ValueError(f"{aligner}: executable not found ({override or ' / '.join(EXECUTABLES[aligner])}). "
                     f"Install it on PATH or pass --{aligner}-exe")


def run_aligner(aligner: str, executable: Path, work: Path, proteins: dict[str, str],
                args: argparse.Namespace) -> tuple[Alignment, dict, str]:
    work.mkdir()
    aliases = tuple(f"seq{i:07d}" for i in range(1, len(proteins) + 1))
    write_fasta(work / "input.fasta", aliases, proteins.values())
    output = work / "aligned.fasta"
    command = [str(executable)]
    mafft_batch = aligner == "mafft" and os.name == "nt" and executable.suffix.lower() in {".bat", ".cmd"}
    if aligner == "mafft":
        command += ["--amino", "--inputorder", "--thread", str(args.threads), *MAFFT_MODES[args.mafft_mode]]
        if mafft_batch:
            command += ["--out", output.name]
        command.append("input.fasta")
    else:
        command += ["-INFILE=input.fasta", "-TYPE=PROTEIN", "-ALIGN", "-OUTPUT=FASTA",
                    f"-OUTFILE={output.name}", "-OUTORDER=INPUT"]

    start = perf_counter()
    result = subprocess.run(command, cwd=work, stdin=subprocess.DEVNULL, capture_output=True, timeout=args.timeout)
    elapsed = perf_counter() - start
    stderr = result.stderr.decode("utf-8", errors="replace")
    stdout_is_alignment = aligner == "mafft" and not mafft_batch
    stdout_log = "" if stdout_is_alignment else result.stdout.decode("utf-8", errors="replace")
    log = "\n".join(part for part in (stdout_log, stderr) if part)
    if result.returncode:
        raise ValueError(f"{aligner} exited with status {result.returncode}\n{log.strip()[-1500:]}")
    if stdout_is_alignment:
        output.write_bytes(result.stdout)
    if not output.is_file():
        raise ValueError(f"{aligner} produced no aligned FASTA\n{log.strip()[-1500:]}")
    try:
        with output.open(encoding="utf-8-sig") as handle:
            aligned = read_fasta(handle)
        by_id = dict(zip(aligned.ids, aligned.sequences))
        if by_id.keys() != set(aliases):
            raise ValueError("output IDs differ from input IDs")
        rows = tuple(by_id[alias].replace(".", "-") for alias in aliases)
        for (name, original), row in zip(proteins.items(), rows):
            if row.replace("-", "") != original:
                raise ValueError(f"output changed the ungapped sequence for {name!r}")
        alignment = Alignment(tuple(proteins), rows)
    except ValueError as exc:
        raise ValueError(f"{aligner}: invalid aligned FASTA ({exc})") from exc

    version = VERSION_PATTERNS[aligner].search(log)
    if version is None and aligner == "mafft":
        try:
            probe = subprocess.run([str(executable), "--version"], cwd=work, stdin=subprocess.DEVNULL, capture_output=True,
                                   timeout=min(args.timeout, 10))
        except (OSError, subprocess.TimeoutExpired):
            pass
        else:
            version = VERSION_PATTERNS[aligner].search((probe.stdout + probe.stderr).decode("utf-8", errors="replace"))
    details = {
        "version": version.group(1) if version else None,
        "command": command,
        "alignment_seconds": elapsed,
        "columns": alignment.length,
        "mean_gap_fraction": sum(row.count("-") for row in rows) / (alignment.size * alignment.length),
    }
    return alignment, details, log


def reference_scores(alignment: Alignment, reference_id: str, score_path: Path) -> list[tuple[str, str, str]]:
    reference = alignment.sequences[alignment.ids.index(reference_id)]
    with score_path.open(encoding="utf-8", newline="") as handle:
        return [(row["position"], row["nongap_coverage"], row["score"])
                for residue, row in zip(reference, csv.DictReader(handle, delimiter="\t"), strict=True)
                if residue not in GAPS]


def write_comparison(path: Path, reference: str, profiles: dict[str, list[tuple[str, str, str]]]) -> dict:
    header = ["reference_position", "reference_residue"]
    for aligner in profiles:
        header.extend(f"{aligner}_{column}" for column in ("column", "nongap_coverage", "score"))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for position, (residue, *values) in enumerate(zip(reference, *profiles.values(), strict=True), 1):
            writer.writerow([position, residue, *(value for site in values for value in site)])
    if len(profiles) == 1:
        return {}
    paired = [(float(mafft[2]), float(clustalw[2]))
              for mafft, clustalw in zip(profiles["mafft"], profiles["clustalw"], strict=True)
              if mafft[2] and clustalw[2]]
    try:
        pearson = correlation(*zip(*paired)) if len(paired) >= 2 else None
    except StatisticsError:
        pearson = None
    return {"jointly_scored_reference_positions": len(paired), "reference_score_pearson": pearson}


def run_comparison(args: argparse.Namespace) -> None:
    data = args.input.read_bytes()
    proteins = read_proteins(data)
    reference_id = args.reference_id or next(iter(proteins))
    if reference_id not in proteins:
        raise ValueError(f"reference ID {reference_id!r} is absent from the input FASTA")
    output = args.output_dir.absolute()
    if output.exists() or output.is_symlink():
        raise ValueError(f"output directory already exists: {output}")
    selected = tuple(EXECUTABLES) if args.aligner == "both" else (args.aligner,)
    executables = {name: resolve_executable(name, getattr(args, f"{name}_exe")) for name in selected}
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="bmi3-msa-", dir=output.parent) as temporary:
        work = Path(temporary)
        results = work / "results"
        results.mkdir()
        score_options = ["--method", args.method, "--max-gap", str(args.max_gap),
                         "--settings-output", str(results / "scoring_settings.json")]
        if args.settings:
            frozen_settings = work / "settings.json"
            frozen_settings.write_bytes(args.settings.read_bytes())
            score_options.extend(("--settings", str(frozen_settings)))
        details = {}
        profiles = {}
        for aligner, executable in executables.items():
            alignment, details[aligner], log = run_aligner(aligner, executable, work / aligner, proteins, args)
            fasta = results / f"{aligner}.fasta"
            scores = results / f"{aligner}.scores.tsv"
            write_fasta(fasta, alignment.ids, alignment.sequences)
            (results / f"{aligner}.log").write_text(log, encoding="utf-8")
            site_options = [*score_options, "--input", str(fasta), "--output", str(scores)]
            if args.method == "scaw-rvne":
                site_options.extend(("--segments-output", str(results / f"{aligner}.segments.tsv")))
            rvne_cli.main(site_options)
            profiles[aligner] = reference_scores(alignment, reference_id, scores)
            details[aligner]["scored_reference_positions"] = sum(bool(site[2]) for site in profiles[aligner])

        comparison = write_comparison(results / "reference_scores.tsv", proteins[reference_id], profiles)
        summary = {"input": str(args.input), "input_sha256": sha256(data).hexdigest(), "sequences": len(proteins),
                   "reference_id": reference_id, "reference_length": len(proteins[reference_id]),
                   "aligners": details, **comparison}
        (results / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        results.rename(output)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="msa_conservation.cli real-data", description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="unaligned protein FASTA")
    parser.add_argument("--output-dir", type=Path, required=True, help="new directory for alignments and scores")
    parser.add_argument("--reference-id", help="first token of the reference FASTA header (default: first sequence)")
    parser.add_argument("--aligner", choices=("both", *EXECUTABLES), default="both")
    parser.add_argument("--mafft-mode", choices=MAFFT_MODES, default="linsi")
    parser.add_argument("--mafft-exe", help="path or executable name (including official Windows mafft.bat)")
    parser.add_argument("--clustalw-exe", help="path or executable name")
    parser.add_argument("--threads", type=int, default=1, help="MAFFT threads (default: 1)")
    parser.add_argument("--timeout", type=int, default=900, help="seconds allowed per aligner (default: 900)")
    parser.add_argument("--method", choices=("rvne", "rvne-optimized", "scaw-rvne"), default="scaw-rvne")
    parser.add_argument("--settings", type=Path, help="frozen RVNE and context parameters")
    parser.add_argument("--max-gap", type=float, default=0.3, help="maximum unweighted gap fraction for scoring")
    args = parser.parse_args(argv)
    if args.threads < 1 or args.timeout < 1 or not 0 <= args.max_gap <= 1:
        parser.error("threads and timeout must be positive and max-gap must be between 0 and 1")
    try:
        run_comparison(args)
    except subprocess.TimeoutExpired as exc:
        parser.exit(2, f"error: {Path(exc.cmd[0]).name} timed out after {args.timeout} seconds\n")
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0

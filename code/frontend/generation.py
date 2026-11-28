"""Web request adapter for the existing INDELible generator (no simulation logic)."""

import inspect
import math
from pathlib import Path
from tempfile import TemporaryDirectory

from data_generation.indelible import Options, find_executable
from data_generation.rvne_simulation import (
    SCENARIOS, SUPPLEMENTAL_SCENARIOS, scenario_rates, simulate, write_simulation,
)


class GenerationError(ValueError):
    pass


BLOCK_SCENARIOS = (
    "blocks", "divergent", "indels", "redundant",
    "conserved-variable-blocks", "three-rate-blocks",
)


def executable_error(exc):
    return str(exc).replace("Set INDELIBLE_BIN or pass --indelible-bin.",
                            "Set INDELIBLE_BIN before starting the local Python service, then restart it.")


def generation_config(max_sequences, max_columns):
    defaults = {key: value.default for key, value in inspect.signature(simulate).parameters.items()}
    engine = Options()
    fields = {
        "size": {"label": "Sequence count", "min": 2, "max": max_sequences, "integer": True},
        "length": {"label": "Ancestor sequence length", "min": 12, "max": max_columns, "integer": True},
        "seed": {"label": "Random seed", "min": 0, "max": 2**31 - 1, "integer": True},
        "divergence": {"label": "Sequence divergence", "min": 0, "max": 3},
        "model": {"label": "Substitution model", "choices": ["WAG", "LG", "JTT"], "default": engine.model},
        "insertion_rate": {"label": "Insertion rate", "min": 0, "max": 1, "default": engine.insertion_rate},
        "deletion_rate": {"label": "Deletion rate", "min": 0, "max": 1, "default": engine.deletion_rate},
        "variable_fraction": {"label": "Variable region fraction", "min": 0, "max": 1,
                              "exclusive": True, "scenarios": BLOCK_SCENARIOS},
        "motif_length": {"label": "Motif length", "min": 1, "max": max_columns // 4,
                         "integer": True, "scenarios": ["motifs", "single-motif"]},
        "redundancy": {"label": "Redundancy fraction", "min": 0, "max": .99,
                       "optional": True, "scenarios": ["redundant"]},
        "gaps": {"label": "Target gap occupancy", "min": 0, "max": .5,
                 "optional": True, "scenarios": ["indels"],
                 "note": "A supplied target calibrates both indel rates, preserving their ratio (1:1 if both are zero)."},
    }
    for key, field in fields.items():
        field.setdefault("default", defaults.get(key))
    try:
        find_executable()
        unavailable = None
    except OSError as exc:
        unavailable = executable_error(exc)
    return {"fields": fields, "scenarios": SCENARIOS + SUPPLEMENTAL_SCENARIOS,
            "scenario": defaults["scenario"], "unavailable": unavailable,
            "maxColumns": max_columns}


def load_generation_options(body, max_sequences, max_columns):
    config = generation_config(max_sequences, max_columns)
    if not isinstance(body, dict):
        raise GenerationError("Provide generation parameters as a JSON object.")
    unknown = set(body) - {"scenario", *config["fields"]}
    if unknown:
        raise GenerationError(f"Unknown generation parameter: {', '.join(sorted(unknown))}.")
    scenario = body.get("scenario", config["scenario"])
    if not isinstance(scenario, str) or scenario not in config["scenarios"]:
        raise GenerationError("Choose a supported simulation scenario.")
    values = {"scenario": scenario}
    for key, field in config["fields"].items():
        applicable = "scenarios" not in field or scenario in field["scenarios"]
        if not applicable and key in body:
            raise GenerationError(f"{field['label']} does not apply to {scenario}.")
        if not applicable:
            continue
        value = body.get(key, field["default"])
        if value is None and field.get("optional"):
            continue
        if "choices" in field:
            valid = isinstance(value, str) and value in field["choices"]
        else:
            valid = (type(value) is int if field.get("integer") else type(value) in (int, float))
            valid = valid and field["min"] <= value <= field["max"] and math.isfinite(value)
            if field.get("exclusive"):
                valid = valid and field["min"] < value < field["max"]
        if not valid:
            choices = ", ".join(field["choices"]) if "choices" in field else (
                f"{'an integer' if field.get('integer') else 'a finite number'} "
                f"{'strictly ' if field.get('exclusive') else ''}between {field['min']} and {field['max']}"
            )
            raise GenerationError(f"{field['label']}: enter {choices}.")
        values[key] = value
    # Reuse the generator's cross-field constraints before starting a job.
    try:
        scenario_rates(scenario, values["length"],
                       values.get("variable_fraction", config["fields"]["variable_fraction"]["default"]),
                       values.get("motif_length", config["fields"]["motif_length"]["default"]))
        options = Options(**{key: values.pop(key) for key in ("model", "insertion_rate", "deletion_rate")})
    except ValueError as exc:
        raise GenerationError(str(exc)) from exc
    return {**values, "options": options}


def generate(options, update, max_columns, max_fasta_bytes):
    update("Locating INDELible")
    try:
        try:
            executable = find_executable()
        except OSError as exc:
            raise GenerationError(executable_error(exc)) from exc
        with TemporaryDirectory(prefix="bmi3-generate-") as directory:
            work = Path(directory)
            update("Simulating protein evolution with INDELible")
            try:
                simulation = simulate(**options, executable=executable, work_dir=work / "run")
            except (OSError, ValueError) as exc:
                # Retain useful engine diagnostics before removing temporary raw files.
                logs = sorted(work.rglob("*.log"), key=lambda path: path.stat().st_mtime, reverse=True)
                details = []
                for path in logs[:4]:
                    tail = path.read_text(encoding="ascii", errors="replace").strip()[-1500:]
                    if tail:
                        details.append(f"{path.parent.name}/{path.name}:\n{tail}")
                raise GenerationError(str(exc) + ("\n" + "\n".join(details) if details else "")) from exc
            if simulation.alignment.length > max_columns:
                raise GenerationError(
                    f"INDELible produced {simulation.alignment.length} aligned columns; this demo accepts "
                    f"up to {max_columns}. Reduce ancestor length, divergence or insertion rate and generate again."
                )
            update("Preparing aligned FASTA and site truth")
            fasta_path, truth_path = work / "aligned.fasta", work / "site_truth.tsv"
            write_simulation(simulation, fasta_path, truth_path, work / "run")
            if fasta_path.stat().st_size > max_fasta_bytes:
                raise GenerationError("Generated FASTA exceeds 1 MiB. Reduce sequence count or ancestor length.")
            fasta = fasta_path.read_text(encoding="ascii")
            truth = truth_path.read_text(encoding="ascii")
            return {"filename": fasta_path.name, "fasta": fasta,
                    "sequenceCount": simulation.alignment.size,
                    "ancestorLength": simulation.metadata["ancestral_length"],
                    "alignedColumns": simulation.alignment.length,
                    "scenario": options["scenario"], "seed": options["seed"],
                    "divergence": simulation.metadata["divergence"]}, {
                        fasta_path.name: fasta, truth_path.name: truth,
                    }
    except (OSError, ValueError) as exc:
        raise GenerationError(f"Generation failed: {exc}") from exc

"""Serve the local conservation workbench and run the existing scoring pipeline."""

import argparse
import csv
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import StringIO
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from threading import Lock, Thread
from time import monotonic
from uuid import uuid4
from zipfile import BadZipFile, ZipFile


FRONTEND = Path(__file__).resolve().parent
CODE = FRONTEND.parent
RVNE = CODE / "rvne"
RVNE_SRC = RVNE / "src"
SCAW_RESULTS = RVNE / "results" / "scaw_revision"
FROZEN = SCAW_RESULTS / "parameters.json"
# The report also archives legacy cohorts and metrics calibrated on legacy families.
BENCHMARK_COHORT = "fresh_indelible"
BENCHMARK_METRICS = {"auprc", "spearman", "boundary_f1", "boundary_f1_exact", "auroc", "scored_fraction"}
MAX_FASTA_BYTES = 1024 * 1024
MAX_REQUEST_BYTES = 2 * MAX_FASTA_BYTES
MAX_SEQUENCES = 500
MAX_COLUMNS = 2000
ALIGNMENT_TIMEOUT = 120
RESULT_TTL = 30 * 60
MAX_RESULTS = 12
STATIC = {
    "/": ("index.html", "text/html"),
    "/styles.css": ("styles.css", "text/css"),
    "/app.js": ("app.js", "text/javascript"),
    "/results.js": ("results.js", "text/javascript"),
    "/benchmarks.js": ("benchmarks.js", "text/javascript"),
    "/ui.js": ("ui.js", "text/javascript"),
    "/generate.js": ("generate.js", "text/javascript"),
}
SYNTHETIC_EXAMPLES = {
    "pattern-c": ("Pattern C – Short conserved motif", 81103),
    "pattern-a": ("Pattern A – Uniform conservation", 81101),
    "pattern-b": ("Pattern B – Central variable block", 81102),
    "pattern-d": ("Pattern D – Three mutation rates", 81105),
    "pattern-d-gaps": ("Pattern D – Three rates with indels", 81106),
}
FASTA_EXAMPLES = {
    "small-alignment": ("Small practice alignment", "aligned", RVNE / "data" / "examples" / "small_alignment.fasta"),
    "globins": ("UniProt globins · Unaligned", "unaligned", RVNE / "data" / "real" / "globins.fasta"),
}

sys.path[:0] = [str(RVNE_SRC), str(CODE)]

from msa_conservation import real_data
from msa_conservation.msa import read_fasta
from msa_conservation.scoring.rvne import DEFAULT_PARAMETERS, RVNEParameters
from msa_conservation.scoring.scaw import ContextParameters, context_from_settings
from generation import GenerationError, generate, generation_config, load_generation_options


def tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def optional_tsv(path):
    return tsv(path) if path.is_file() else []


def number(value):
    if value in (None, ""):
        return None
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def source(path):
    return path.relative_to(CODE.parent).as_posix()


def example_labels():
    try:
        with ZipFile(SCAW_RESULTS / "fresh_inputs.zip") as archive:
            members = set(archive.namelist())
    except (OSError, BadZipFile):
        members = set()
    available = {key: label for key, (label, seed) in SYNTHETIC_EXAMPLES.items()
                 if {f"{seed}/alignment.fasta", f"{seed}/truth.tsv"} <= members}
    return {**available, **{key: entry[0] for key, entry in FASTA_EXAMPLES.items()}}


def example(key):
    if key in SYNTHETIC_EXAMPLES:
        label, seed = SYNTHETIC_EXAMPLES[key]
        with ZipFile(SCAW_RESULTS / "fresh_inputs.zip") as archive:
            fasta = archive.read(f"{seed}/alignment.fasta").decode("ascii")
            rows = csv.DictReader(StringIO(archive.read(f"{seed}/truth.tsv").decode("ascii")), delimiter="\t")
            truth = [{"position": int(row["position"]), "rate": float(row["rate"]),
                      "conserved": row["conserved"] == "1"} for row in rows]
        return {"filename": f"{key}-{seed}.fasta", "mode": "aligned", "fasta": fasta, "truth": truth,
                "source": f"{label} – INDELible 1.03, seed {seed}"}
    label, mode, path = FASTA_EXAMPLES[key]
    return {"filename": path.name, "mode": mode, "fasta": path.read_text(encoding="utf-8"),
            "truth": [], "source": f"{label} · {source(path)}"}


def globin_case():
    directory = RVNE / "results" / "globins"
    path = directory / "reference_scores.tsv"
    summary_path = directory / "summary.json"
    if not path.is_file() or not summary_path.is_file():
        return None
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("reference_id") != "P02144" or not {"mafft", "clustalw"} <= summary["aligners"].keys():
        return None
    sites = [{"position": int(row["reference_position"]), "mafft": number(row["mafft_score"]),
              "clustalw": number(row["clustalw_score"])} for row in tsv(path)]
    if [site["position"] for site in sites] != list(range(1, summary["reference_length"] + 1)):
        raise ValueError("Globin scores do not cover each P02144 reference position in order.")
    return {"sites": sites, "source": source(path), "pearson": summary.get("reference_score_pearson"),
            "joint": sum(site["mafft"] is not None and site["clustalw"] is not None for site in sites)}


def benchmarks():
    summary = SCAW_RESULTS / "summary.tsv"
    baselines = SCAW_RESULTS / "baseline_summary.tsv"
    runtime = SCAW_RESULTS / "runtime.tsv"
    datasets = {}
    for row in optional_tsv(summary):
        if row["cohort"] != BENCHMARK_COHORT or row["metric"] not in BENCHMARK_METRICS:
            continue
        dataset = datasets.setdefault(row["cohort"], {"rows": [], "groups": {}})
        dataset["rows"].append({
            "method": row["method"], "pattern": row["scenario"], "metric": row["metric"],
            "mean": number(row["mean"]), "sd": number(row["sd"]), "n": int(row["n"]),
        })
        dataset["groups"].setdefault(row["scenario"], set()).add(row["method"])
    for dataset in datasets.values():
        dataset["groups"] = {key: sorted(methods) for key, methods in dataset["groups"].items()}
    comparison_rows = [
        {"method": row["method"], "metric": row["metric"], "n": int(row["n"]),
         "mean": number(row["mean"]), "sd": number(row["sd"])}
        for row in optional_tsv(baselines)
        if row["cohort"] == BENCHMARK_COHORT and row["scenario"] == "all" and row["metric"] in BENCHMARK_METRICS
    ]
    try:
        case = globin_case()
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"Globin results could not be loaded: {exc}", file=sys.stderr)
        case = {"error": "The recorded globin outputs could not be read. Regenerate the MAFFT and ClustalW comparison."}
    return {
        "datasets": datasets,
        "source": source(summary),
        "baselines": comparison_rows,
        "baseline_source": source(baselines),
        "runtime": [{"method": row["method"], "columns": int(row["columns"]),
                     "seconds": number(row["median_seconds"]), "repeats": int(row["repeats"])}
                    for row in optional_tsv(runtime) if row["median_seconds"]],
        "runtime_source": source(runtime),
        "case": case,
    }


class AnalysisError(ValueError):
    pass


def aligners():
    available = {}
    for name in real_data.EXECUTABLES:
        try:
            available[name] = real_data.resolve_executable(name, os.environ.get(f"BMI3_{name.upper()}"))
        except ValueError:
            available[name] = None
    return available


def custom_settings(body, method):
    if not isinstance(body, dict) or not isinstance(rvne := body.get("rvne"), dict):
        raise AnalysisError("Provide custom RVNE parameters.")
    rvne_values = {key: rvne.get(key) for key in vars(DEFAULT_PARAMETERS)}
    if any(type(value) not in (int, float) for value in rvne_values.values()):
        raise AnalysisError("Enter numeric values for r, s1 and s2.")
    settings = {"rvne": vars(RVNEParameters(**rvne_values))}
    if method != "scaw-rvne":
        return settings

    context = body.get("context")
    if not isinstance(context, dict):
        raise AnalysisError("Provide custom SCAW parameters.")
    if set(context) != {"penalty_scale"}:
        raise AnalysisError("SCAW v2 requires only penalty_scale. Reload the page to update the settings form.")
    settings["context"] = vars(context_from_settings(context))
    return settings


def load_options(body):
    if not isinstance(body, dict):
        raise AnalysisError("Provide an analysis request.")
    mode = body.get("mode")
    method = body.get("method")
    preset = body.get("preset")
    aligner = body.get("aligner", "mafft")
    if mode not in ("aligned", "unaligned") or method not in ("scaw-rvne", "rvne"):
        raise AnalysisError("Choose a valid input mode and scoring method.")
    if preset not in ("frozen", "default", "custom") or aligner not in real_data.EXECUTABLES:
        raise AnalysisError("Choose valid scoring settings and alignment tool.")
    gap = body.get("maxGap")
    if type(gap) not in (int, float) or not 0 <= gap <= 1:
        raise AnalysisError("Maximum gap fraction must be between 0 and 1.")
    key = body.get("example")
    if key is not None:
        if not isinstance(key, str) or key not in example_labels():
            raise AnalysisError("Unknown example.")
        chosen = example(key)
        if mode != chosen["mode"]:
            raise AnalysisError("Use the input mode indicated for this example.")
    else:
        filename = body.get("filename")
        fasta = body.get("fasta")
        if not isinstance(filename, str) or not isinstance(fasta, str) or not fasta.strip():
            raise AnalysisError("Upload a nonempty protein FASTA.")
        if len(fasta.encode("utf-8")) > MAX_FASTA_BYTES:
            raise AnalysisError("This local demo accepts FASTA files up to 1 MiB.")
        chosen = {"filename": filename.replace("\\", "/").rsplit("/", 1)[-1][:180] or "alignment.fasta",
                  "fasta": fasta, "truth": [], "source": None}
    try:
        settings = custom_settings(body.get("parameters"), method) if preset == "custom" else None
    except ValueError as exc:
        raise AnalysisError(str(exc)) from exc
    return {**chosen, "example": key, "mode": mode, "method": method, "preset": preset, "aligner": aligner,
            "max_gap": gap, "parameters": settings}


def check_dimensions(size, length):
    if size < 2:
        raise AnalysisError("FASTA format: provide at least two sequences with unique IDs.")
    if size > MAX_SEQUENCES or length > MAX_COLUMNS:
        raise AnalysisError(f"This demo supports up to {MAX_SEQUENCES} sequences and {MAX_COLUMNS} columns.")


def run_scoring(alignment, options, work):
    fasta = work / "alignment.fasta"
    scores = work / "scores.tsv"
    segments = work / "segments.tsv"
    settings = work / "settings.json"
    real_data.write_fasta(fasta, alignment.ids, alignment.sequences)
    command = [
        sys.executable, "-m", "msa_conservation.cli", "score", "--input", str(fasta),
        "--method", options["method"], "--max-gap", str(options["max_gap"]),
        "--output", str(scores), "--settings-output", str(settings),
    ]
    if options["preset"] == "frozen":
        command.extend(("--settings", str(FROZEN)))
    elif options["preset"] == "custom":
        custom = work / "custom-settings.json"
        custom.write_text(json.dumps(options["parameters"]), encoding="utf-8")
        command.extend(("--settings", str(custom)))
    if options["method"] == "scaw-rvne":
        command.extend(("--segments-output", str(segments)))
    try:
        process = subprocess.run(
            command, cwd=RVNE, env={**os.environ, "PYTHONPATH": str(RVNE_SRC), "PYTHONIOENCODING": "utf-8"},
            stdin=subprocess.DEVNULL, capture_output=True, encoding="utf-8", errors="replace", timeout=120,
        )
    except subprocess.TimeoutExpired as exc:
        raise AnalysisError("Scoring timed out. Try a smaller alignment.") from exc
    if process.returncode:
        details = process.stderr.strip().splitlines()
        raise AnalysisError("Scoring failed. " + (details[-1] if details else "Check the scoring environment."))
    sites = [{key: number(value) for key, value in row.items()} for row in tsv(scores)]
    fitted = [{**{key: int(row[key]) for key in ("segment_id", "start", "end", "valid_columns")},
               "raw_mean": number(row["raw_mean"]), "class": row["class"]} for row in optional_tsv(segments)]
    downloads = {
        "alignment.fasta": fasta.read_text(encoding="utf-8"),
        "scores.tsv": scores.read_text(encoding="utf-8"),
        "settings.json": settings.read_text(encoding="utf-8"),
    }
    if segments.is_file():
        downloads["segments.tsv"] = segments.read_text(encoding="utf-8")
    return sites, fitted, downloads


def analyze(options, update):
    aligned = options["mode"] == "aligned"
    update("Reading FASTA")
    try:
        if aligned:
            alignment = read_fasta(options["fasta"].removeprefix("\ufeff").splitlines())
        else:
            proteins = real_data.read_proteins(options["fasta"].encode("utf-8"))
    except (ValueError, UnicodeError) as exc:
        raise AnalysisError(f"FASTA format: {exc}") from exc
    aligner = options["aligner"] if not aligned else None
    aligner_details = None
    with TemporaryDirectory(prefix="bmi3-demo-") as directory:
        work = Path(directory)
        if aligner:
            check_dimensions(len(proteins), max(map(len, proteins.values())))
            executable = aligners()[aligner]
            if executable is None:
                label = {"mafft": "MAFFT", "clustalw": "ClustalW"}[aligner]
                raise AnalysisError(f"{label} is unavailable. Install it on PATH or use an already aligned FASTA.")
            update("Aligning")
            config = argparse.Namespace(threads=1, mafft_mode="auto", timeout=ALIGNMENT_TIMEOUT)
            try:
                alignment, aligner_details, _ = real_data.run_aligner(aligner, executable, work / aligner, proteins, config)
            except subprocess.TimeoutExpired as exc:
                raise AnalysisError("Alignment timed out. Try a smaller family.") from exc
            except (OSError, ValueError) as exc:
                raise AnalysisError(f"Alignment failed: {exc}") from exc
        check_dimensions(alignment.size, alignment.length)
        update("Scoring")
        sites, segments, downloads = run_scoring(alignment, options, work)
    reference = "P02144" if "P02144" in alignment.ids and options["source"] else None
    reference_positions = []
    if reference:
        position = 0
        for residue in alignment.sequences[alignment.ids.index(reference)]:
            if residue not in "-.":
                position += 1
                reference_positions.append(position)
            else:
                reference_positions.append(None)
    return {
        "filename": options["filename"], "method": options["method"], "preset": options["preset"],
        "example": options["example"],
        "aligner": aligner, "alignerVersion": aligner_details["version"] if aligner_details else None,
        "ids": alignment.ids, "sequences": alignment.sequences, "sites": sites, "segments": segments,
        "scored": sum(site["score"] is not None for site in sites), "source": options["source"],
        "truth": options["truth"], "reference": reference, "referencePositions": reference_positions,
    }, downloads


class Workbench(ThreadingHTTPServer):
    def __init__(self, port):
        super().__init__(("127.0.0.1", port), Handler)
        self.jobs = {}
        self.lock = Lock()

    def cleanup(self):
        expired = [key for key, job in self.jobs.items()
                   if job.get("finished", float("inf")) < monotonic() - RESULT_TTL]
        for key in expired:
            del self.jobs[key]

    def submit(self, options, kind="analysis"):
        with self.lock:
            self.cleanup()
            if any("finished" not in job for job in self.jobs.values()):
                return None
            if len(self.jobs) >= MAX_RESULTS:
                del self.jobs[next(iter(self.jobs))]
            key = uuid4().hex
            self.jobs[key] = {"state": "running", "stage": "Preparing generation" if kind == "generation" else "Reading FASTA"}
        Thread(target=self.run_job, args=(key, options, kind), daemon=True).start()
        return key

    def run_job(self, key, options, kind="analysis"):
        def update(stage):
            with self.lock:
                self.jobs[key]["stage"] = stage

        try:
            result, downloads = (generate(options, update, MAX_COLUMNS, MAX_FASTA_BYTES)
                                 if kind == "generation" else analyze(options, update))
            job = {"state": "done", "result": result, "downloads": downloads}
        except (AnalysisError, GenerationError) as exc:
            job = {"state": "error", "error": str(exc)}
        except Exception as exc:
            print(f"{kind.capitalize()} error: {exc}", file=sys.stderr)
            message = f"Generation failed: {exc}" if kind == "generation" else "Scoring failed. Check the server terminal for details."
            job = {"state": "error", "error": message}
        with self.lock:
            self.jobs[key] = {**job, "finished": monotonic()}


class Handler(BaseHTTPRequestHandler):
    def send_data(self, status, content, kind, download=None):
        self.send_response(status)
        self.send_header("Content-Type", f"{kind}; charset=utf-8")
        self.send_header("Content-Length", len(content))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy",
                         "default-src 'self'; style-src 'self' 'unsafe-inline'; frame-ancestors 'none'")
        if download:
            self.send_header("Content-Disposition", f'attachment; filename="{download}"')
        self.end_headers()
        self.wfile.write(content)

    def send_json(self, value, status=200):
        payload = json.dumps(value, ensure_ascii=True, allow_nan=False).encode("utf-8")
        self.send_data(status, payload, "application/json")

    def local_request(self):
        port = self.server.server_port
        hosts = {f"127.0.0.1:{port}", f"localhost:{port}"}
        if port == 80:
            hosts.update(("127.0.0.1", "localhost"))
        origin = self.headers.get("Origin")
        valid = self.headers.get("Host") in hosts and (origin is None or origin in {f"http://{h}" for h in hosts})
        if not valid:
            self.send_json({"error": "Open the local workbench URL printed in the terminal."}, 403)
        return valid

    def do_GET(self):
        if not self.local_request():
            return
        try:
            self.route_get()
        except (OSError, ValueError, KeyError) as exc:
            print(f"Read error: {exc}", file=sys.stderr)
            self.send_json({"error": "The requested project data could not be loaded."}, 500)

    def route_get(self):
        if self.path in STATIC:
            filename, kind = STATIC[self.path]
            self.send_data(200, (FRONTEND / filename).read_bytes(), kind)
        elif self.path == "/api/config":
            frozen = json.loads(FROZEN.read_text(encoding="utf-8")) if FROZEN.is_file() else None
            self.send_json({"aligners": {name: path is not None for name, path in aligners().items()},
                            "frozen": {key: frozen[key] for key in ("rvne", "context")} if frozen else None,
                            "defaults": {"rvne": vars(DEFAULT_PARAMETERS), "context": vars(ContextParameters())},
                            "examples": example_labels()})
        elif self.path.startswith("/api/examples/"):
            key = self.path.removeprefix("/api/examples/")
            if key in example_labels():
                self.send_json({name: value for name, value in example(key).items() if name != "truth"})
            else:
                self.send_json({"error": "Unknown example."}, 404)
        elif self.path == "/api/benchmarks":
            self.send_json(benchmarks())
        elif self.path == "/api/generation/config":
            self.send_json(generation_config(MAX_SEQUENCES, MAX_COLUMNS))
        elif self.path.startswith("/api/jobs/"):
            self.get_job(self.path.removeprefix("/api/jobs/").split("/"))
        else:
            self.send_json({"error": "Page not found."}, 404)

    def get_job(self, parts):
        with self.server.lock:
            self.server.cleanup()
            job = self.server.jobs.get(parts[0])
        if job is None:
            self.send_json({"error": "This result expired. Run it again."}, 404)
        elif len(parts) == 1:
            self.send_json({key: value for key, value in job.items() if key not in {"downloads", "finished"}})
        elif len(parts) == 3 and parts[1] == "download" and parts[2] in job.get("downloads", {}):
            filename = parts[2]
            kind = {"tsv": "text/tab-separated-values", "json": "application/json"}.get(filename.split(".")[-1], "text/plain")
            self.send_data(200, job["downloads"][filename].encode("utf-8"), kind, download=filename)
        else:
            self.send_json({"error": "Download not available."}, 404)

    def do_POST(self):
        if not self.local_request():
            return
        if self.path not in ("/api/jobs", "/api/generate"):
            self.send_json({"error": "Page not found."}, 404)
            return
        if self.headers.get_content_type() != "application/json":
            self.send_json({"error": "Send a JSON request."}, 415)
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= MAX_REQUEST_BYTES:
                raise AnalysisError("Upload a protein FASTA of up to 1 MiB.")
            body = json.loads(self.rfile.read(length))
            options = (load_generation_options(body, MAX_SEQUENCES, MAX_COLUMNS)
                       if self.path == "/api/generate" else load_options(body))
        except (AnalysisError, GenerationError) as exc:
            self.send_json({"error": str(exc)}, 400)
            return
        except (ValueError, UnicodeError, TypeError):
            self.send_json({"error": "Invalid JSON request."}, 400)
            return
        key = self.server.submit(options, "generation" if self.path == "/api/generate" else "analysis")
        if key is None:
            self.send_json({"error": "An analysis or generation is already running. Try again once it finishes."}, 409)
        else:
            self.send_json({"id": key}, 202)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    with Workbench(args.port) as server:
        print(f"Protein Conservation Scorer: http://127.0.0.1:{server.server_port}/", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()

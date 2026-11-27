"""Serve the local conservation workbench and run the existing scoring pipeline."""

import argparse
import csv
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
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


FRONTEND = Path(__file__).resolve().parent
CODE = FRONTEND.parent
RVNE = CODE / "rvne"
RVNE_SRC = RVNE / "src"
REFERENCE = RVNE / "results" / "reference_run"
SUPPLEMENTAL = RVNE / "results" / "supplemental_run"
FROZEN = REFERENCE / "parameters.json"
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
}
SYNTHETIC_EXAMPLES = {
    "pattern-c": ("Pattern C · Short conserved motif", "c_50124.fasta",
                  {"scenario": "single-motif", "seed": 50124, "motif_length": 3}),
    "pattern-a": ("Pattern A · Uniform conservation", "a_50104.fasta",
                  {"scenario": "conserved", "seed": 50104}),
    "pattern-b": ("Pattern B · Central variable block", "b_50116.fasta",
                  {"scenario": "conserved-variable-blocks", "seed": 50116}),
    "pattern-d": ("Pattern D · Three mutation rates", "d_50140.fasta",
                  {"scenario": "three-rate-blocks", "seed": 50140}),
    "pattern-d-gaps": ("Pattern D · Three rates with gaps", "d_50152.fasta",
                       {"scenario": "three-rate-blocks", "seed": 50140,
                        "gap_by_rate": (0.03, 0.12, 0.25), "gap_seed": 50152}),
}
FASTA_EXAMPLES = {
    "small-alignment": ("Small practice alignment", "aligned", RVNE / "data" / "examples" / "small_alignment.fasta"),
    "globins": ("UniProt globins · Unaligned", "unaligned", RVNE / "data" / "real" / "globins.fasta"),
}
EXAMPLES = {key: entry[0] for key, entry in (SYNTHETIC_EXAMPLES | FASTA_EXAMPLES).items()}

sys.path[:0] = [str(RVNE_SRC), str(CODE)]

from msa_conservation import real_data
from msa_conservation.msa import read_fasta
from msa_conservation.scoring.rvne import DEFAULT_PARAMETERS, RVNEParameters
from msa_conservation.scoring.scaw import ContextParameters
from data_generation.legacy_reference import simulate


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


@lru_cache(maxsize=len(EXAMPLES))
def example(key):
    if key in SYNTHETIC_EXAMPLES:
        label, filename, settings = SYNTHETIC_EXAMPLES[key]
        sampled = simulate(**settings)
        fasta = "".join(f">{name}\n{sequence}\n" for name, sequence in
                        zip(sampled.alignment.ids, sampled.alignment.sequences))
        truth = [{"position": i, "rate": float(rate), "conserved": bool(conserved)}
                 for i, (rate, conserved) in enumerate(zip(sampled.rates, sampled.conserved), 1)]
        return {"filename": filename, "mode": "aligned", "fasta": fasta, "truth": truth,
                "source": f"{label} · seed {settings.get('gap_seed', settings['seed'])}"}
    label, mode, path = FASTA_EXAMPLES[key]
    return {"filename": path.name, "mode": mode, "fasta": path.read_text(encoding="utf-8"),
            "truth": [], "source": f"{label} · {source(path)}"}


def benchmark_dataset(summary, families, metric_aliases=None):
    if not summary.is_file():
        return {"rows": [], "available": False, "source": source(summary)}
    aliases = metric_aliases or {}
    observations = [
        {"method": row["method"], "pattern": row["scenario"], "metric": aliases.get(row["metric"], row["metric"]),
         "mean": number(row["mean"]), "sd": number(row["sd"]), "n": int(row["n"])}
        for row in tsv(summary)
    ]
    groups = {}
    for row in optional_tsv(families):
        groups.setdefault(row["scenario"], set()).add(row["method"])
    for row in observations:
        if row["pattern"] != "all":
            groups.setdefault(row["pattern"], set()).add(row["method"])
    return {"rows": observations, "available": True, "source": source(summary),
            "families_source": source(families), "groups": {key: sorted(methods) for key, methods in groups.items()}}


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
    return {"sites": sites, "source": source(path),
            "joint": sum(site["mafft"] is not None and site["clustalw"] is not None for site in sites)}


def benchmarks():
    settings = SUPPLEMENTAL / "settings_summary.tsv"
    runtime = REFERENCE / "runtime.tsv"
    return {
        "datasets": {
            "supplemental": benchmark_dataset(SUPPLEMENTAL / "summary.tsv", SUPPLEMENTAL / "families.tsv"),
            "reference": benchmark_dataset(
                REFERENCE / "test_summary.tsv", REFERENCE / "test_families.tsv", {"specificity": "variable_recall"},
            ),
        },
        "motifs": [{"method": row["method"], "width": int(row["scenario"].split("=")[1]),
                    "mean": number(row["mean"]), "n": int(row["n"])}
                   for row in optional_tsv(settings)
                   if row["scenario"].startswith("C motif_length=") and row["metric"] == "sensitivity"
                   and row["method"] in {"E2", "E3", "E7", "E8", "Shannon"}],
        "motif_source": source(settings),
        "runtime": [{"stage": row["stage"], "sequences": int(row["N"]), "columns": int(row["L"]),
                     "seconds": number(row["median_seconds"]), "repeats": int(row["repeats"])}
                    for row in optional_tsv(runtime) if row["median_seconds"]],
        "runtime_source": source(runtime),
        "case": globin_case(),
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
    weights = {key: context.get(key) for key in ("beta", "smoothing", "min_support")}
    if any(type(value) not in (int, float) or not 0 <= value <= sys.float_info.max
           for value in weights.values()):
        raise AnalysisError("Segmentation penalty, window smoothing and minimum support must be finite and nonnegative.")
    radii = context.get("radii")
    if not isinstance(radii, list) or any(type(radius) is not int or not 0 <= radius <= MAX_COLUMNS for radius in radii):
        raise AnalysisError(f"Window radii must be whole numbers from 0 to {MAX_COLUMNS}.")
    minimum = context.get("min_segment")
    maximum = context.get("max_segment")
    if type(minimum) is not int or (maximum is not None and type(maximum) is not int):
        raise AnalysisError("Enter a whole number for each segment length or leave the maximum blank.")
    settings["context"] = vars(ContextParameters(**weights, radii=tuple(radii), min_segment=minimum,
                                                  max_segment=maximum))
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
        if not isinstance(key, str) or key not in EXAMPLES:
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

    def submit(self, options):
        with self.lock:
            self.cleanup()
            if any("finished" not in job for job in self.jobs.values()):
                return None
            if len(self.jobs) >= MAX_RESULTS:
                del self.jobs[next(iter(self.jobs))]
            key = uuid4().hex
            self.jobs[key] = {"state": "running", "stage": "Reading FASTA"}
        Thread(target=self.run_job, args=(key, options), daemon=True).start()
        return key

    def run_job(self, key, options):
        def update(stage):
            with self.lock:
                self.jobs[key]["stage"] = stage

        try:
            result, downloads = analyze(options, update)
            job = {"state": "done", "result": result, "downloads": downloads}
        except AnalysisError as exc:
            job = {"state": "error", "error": str(exc)}
        except Exception as exc:
            print(f"Analysis error: {exc}", file=sys.stderr)
            job = {"state": "error", "error": "Scoring failed. Check the server terminal for details."}
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
                            "examples": EXAMPLES})
        elif self.path.startswith("/api/examples/"):
            key = self.path.removeprefix("/api/examples/")
            if key in EXAMPLES:
                self.send_json({name: value for name, value in example(key).items() if name != "truth"})
            else:
                self.send_json({"error": "Unknown example."}, 404)
        elif self.path == "/api/benchmarks":
            self.send_json(benchmarks())
        elif self.path.startswith("/api/jobs/"):
            self.get_job(self.path.removeprefix("/api/jobs/").split("/"))
        else:
            self.send_json({"error": "Page not found."}, 404)

    def get_job(self, parts):
        with self.server.lock:
            self.server.cleanup()
            job = self.server.jobs.get(parts[0])
        if job is None:
            self.send_json({"error": "This analysis expired. Run it again."}, 404)
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
        if self.path != "/api/jobs":
            self.send_json({"error": "Page not found."}, 404)
            return
        if self.headers.get_content_type() != "application/json":
            self.send_json({"error": "Upload a JSON analysis request."}, 415)
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= MAX_REQUEST_BYTES:
                raise AnalysisError("Upload a protein FASTA of up to 1 MiB.")
            options = load_options(json.loads(self.rfile.read(length)))
        except AnalysisError as exc:
            self.send_json({"error": str(exc)}, 400)
            return
        except (ValueError, UnicodeError, TypeError):
            self.send_json({"error": "Invalid analysis request."}, 400)
            return
        key = self.server.submit(options)
        if key is None:
            self.send_json({"error": "An analysis is already running. Try again once it finishes."}, 409)
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

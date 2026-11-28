"""Run with python -m unittest discover -s code/frontend/tests -v."""

import csv
from io import StringIO
import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from threading import Event, Thread
from time import monotonic, sleep
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from urllib.error import HTTPError
from urllib.request import Request, urlopen
from zipfile import ZipFile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import server
from generation import GenerationError, generate, generation_config, load_generation_options
from data_generation.indelible import Options, find_executable
from data_generation.rvne_simulation import SCENARIOS, SUPPLEMENTAL_SCENARIOS, simulate
from msa_conservation.msa import read_fasta


class GenerationValidationTests(unittest.TestCase):
    def options(self, body):
        return load_generation_options(body, server.MAX_SEQUENCES, server.MAX_COLUMNS)

    def test_example_reads_precomputed_indelible_alignment_and_truth(self):
        with TemporaryDirectory() as directory:
            folder = Path(directory)
            with patch.object(server, "SCAW_RESULTS", folder):
                self.assertNotIn("pattern-c", server.example_labels())
                with ZipFile(folder / "fresh_inputs.zip", "w") as archive:
                    archive.writestr("81103/alignment.fasta", ">a\nAS\n>b\nAT\n")
                    archive.writestr("81103/truth.tsv", "position\trate\tconserved\n1\t0.06\t1\n2\t1.0\t0\n")
                self.assertIn("pattern-c", server.example_labels())
                example = server.example("pattern-c")
                self.assertIn("INDELible 1.03, seed 81103", example["source"])
                self.assertEqual(example["fasta"], ">a\nAS\n>b\nAT\n")
                self.assertEqual(example["truth"], [{"position": 1, "rate": .06, "conserved": True},
                                                    {"position": 2, "rate": 1., "conserved": False}])

    def test_defaults_follow_generator(self):
        import inspect
        options = self.options({})
        defaults = inspect.signature(simulate).parameters
        for key in ("size", "length", "seed", "scenario", "divergence", "variable_fraction"):
            self.assertEqual(options[key], defaults[key].default)
        self.assertEqual(options["options"], Options())
        config = generation_config(500, 2000)
        self.assertEqual(tuple(config["scenarios"]), SCENARIOS + SUPPLEMENTAL_SCENARIOS)

    def test_invalid_and_inapplicable_parameters(self):
        for body in (None, [], {"size": True}, {"size": 1}, {"size": 501}, {"size": 2.5},
                     {"length": 11}, {"seed": -1}, {"seed": 10**1000},
                     {"divergence": float("nan")}, {"divergence": float("inf")},
                     {"model": "DNA"}, {"model": []}, {"insertion_rate": -1},
                     {"deletion_rate": 1.01}, {"scenario": []}, {"scenario": "custom"},
                     {"variable_fraction": 1}, {"scenario": "conserved", "motif_length": 3},
                     {"scenario": "single-motif", "length": 12, "motif_length": 4},
                     {"scenario": "redundant", "redundancy": 1},
                     {"scenario": "indels", "gaps": .6}, {"executable": "arbitrary-path"}):
            with self.subTest(body=repr(body)[:120]):
                with self.assertRaises(GenerationError):
                    self.options(body)

    def test_all_scenarios_and_conditional_options(self):
        for scenario in SCENARIOS + SUPPLEMENTAL_SCENARIOS:
            with self.subTest(scenario=scenario):
                self.assertEqual(self.options({"scenario": scenario})["scenario"], scenario)
        self.assertEqual(self.options({"scenario": "single-motif", "motif_length": 3})["motif_length"], 3)
        self.assertEqual(self.options({"scenario": "indels", "gaps": .1})["gaps"], .1)
        self.assertNotIn("redundancy", self.options({"scenario": "redundant"}))

    def test_generated_width_and_file_size_limits(self):
        options = self.options({})
        too_wide = SimpleNamespace(alignment=SimpleNamespace(length=2001))
        with patch("generation.find_executable", return_value=Path("test.exe")):
            with patch("generation.simulate", autospec=True, return_value=too_wide):
                with self.assertRaisesRegex(GenerationError, "2001 aligned columns"):
                    generate(options, lambda stage: None, 2000, 1024)

            def write_large(simulation, fasta, truth, work):
                fasta.write_text("A" * 1025, encoding="ascii")

            with patch("generation.simulate", autospec=True, return_value=SimpleNamespace(
                    alignment=SimpleNamespace(length=240))):
                with patch("generation.write_simulation", side_effect=write_large):
                    with self.assertRaisesRegex(GenerationError, "exceeds 1 MiB"):
                        generate(options, lambda stage: None, 2000, 1024)

    def test_simulator_diagnostics_survive_temporary_cleanup(self):
        def failed_simulation(*args, work_dir, **kwargs):
            work_dir.mkdir()
            (work_dir / "stderr.log").write_text("specific simulator error", encoding="ascii")
            raise ValueError("INDELible exited with status 2")

        options = self.options({})
        with patch("generation.find_executable", return_value=Path("test.exe")):
            with patch("generation.simulate", autospec=True, side_effect=failed_simulation):
                with self.assertRaisesRegex(GenerationError, "specific simulator error"):
                    generate(options, lambda stage: None, 2000, 1024)


class GenerationHTTPTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.service = server.Workbench(0)
        cls.thread = Thread(target=cls.service.serve_forever, daemon=True)
        cls.thread.start()
        cls.url = f"http://127.0.0.1:{cls.service.server_port}"

    @classmethod
    def tearDownClass(cls):
        cls.service.shutdown()
        cls.service.server_close()
        cls.thread.join()

    def request(self, path, body=None):
        request = Request(self.url + path, data=json.dumps(body).encode() if body is not None else None,
                          headers={"Content-Type": "application/json"})
        try:
            with urlopen(request, timeout=10) as response:
                return response.status, response.read(), response.headers
        except HTTPError as exc:
            return exc.code, exc.read(), exc.headers

    def job(self, body, path="/api/generate"):
        status, raw, _ = self.request(path, body)
        self.assertEqual(status, 202, raw)
        key = json.loads(raw)["id"]
        deadline = monotonic() + 90
        while monotonic() < deadline:
            status, raw, _ = self.request(f"/api/jobs/{key}")
            self.assertEqual(status, 200)
            job = json.loads(raw)
            if job["state"] != "running":
                return key, job
            sleep(.1)
        self.fail("Job did not finish within 90 seconds")

    def test_static_config_and_invalid_requests(self):
        for path in ("/", "/app.js", "/generate.js", "/styles.css", "/api/config", "/api/generation/config"):
            with self.subTest(path=path):
                self.assertEqual(self.request(path)[0], 200)
        for body in ({"size": True}, {"scenario": "single-motif", "length": 12, "motif_length": 4}):
            status, raw, _ = self.request("/api/generate", body)
            self.assertEqual(status, 400)
            self.assertTrue(json.loads(raw)["error"])

    def test_missing_indelible_config_and_job(self):
        with patch.dict(os.environ, {"INDELIBLE_BIN": "bmi3-deliberately-missing-indelible.exe"}):
            config = json.loads(self.request("/api/generation/config")[1])
            self.assertIn("INDELIBLE_BIN", config["unavailable"])
            _, job = self.job({})
            self.assertEqual(job["state"], "error")
            self.assertIn("executable not found", job["error"])
            self.assertIn("INDELIBLE_BIN", job["error"])

    def test_simulator_failure(self):
        with patch("generation.simulate", autospec=True, side_effect=ValueError("INDELible timed out after 60s")):
            with patch("generation.find_executable", return_value=Path("test.exe")):
                _, job = self.job({})
        self.assertEqual(job["state"], "error")
        self.assertIn("timed out after 60s", job["error"])

    def test_busy_job_is_not_evicted_and_download_expires(self):
        started, release = Event(), Event()

        def paused(*args):
            started.set()
            release.wait(10)
            raise GenerationError("deliberately stopped test job")

        with patch("server.generate", side_effect=paused):
            status, raw, _ = self.request("/api/generate", {})
            self.assertEqual(status, 202)
            key = json.loads(raw)["id"]
            try:
                self.assertTrue(started.wait(2))
                self.assertEqual(self.request("/api/generate", {})[0], 409)
                self.assertEqual(self.request(f"/api/jobs/{key}/download/aligned.fasta")[0], 404)
            finally:
                release.set()
            deadline = monotonic() + 3
            while monotonic() < deadline:
                with self.service.lock:
                    if "finished" in self.service.jobs[key]:
                        break
                sleep(.01)
        with self.service.lock:
            self.service.jobs[key]["finished"] = monotonic() - server.RESULT_TTL - 1
        self.assertEqual(self.request(f"/api/jobs/{key}")[0], 404)

    def test_real_default_and_modified_generation_score_and_truth(self):
        try:
            find_executable()
        except FileNotFoundError:
            self.skipTest("INDELible 1.03 is required for real generation")
        cases = [({}, (50, 240)),
                 ({"scenario": "three-rate-blocks", "size": 8, "length": 80, "seed": 43,
                   "divergence": 1.1, "model": "LG", "insertion_rate": .08,
                   "deletion_rate": .05, "variable_fraction": .25}, (8, 80)),
                 ({"scenario": "single-motif", "size": 4, "length": 60, "seed": 9,
                   "motif_length": 3, "model": "JTT"}, (4, 60))]
        for parameters, (size, ancestor) in cases:
            with self.subTest(parameters=parameters):
                key, job = self.job(parameters)
                self.assertEqual(job["state"], "done", job)
                result = job["result"]
                alignment = read_fasta(result["fasta"].splitlines())
                self.assertEqual((alignment.size, result["ancestorLength"]), (size, ancestor))
                self.assertEqual(alignment.length, result["alignedColumns"])
                code, data, headers = self.request(f"/api/jobs/{key}/download/aligned.fasta")
                self.assertEqual(code, 200)
                self.assertIn("attachment", headers["Content-Disposition"])
                self.assertEqual(data.decode(), result["fasta"])
                code, data, _ = self.request(f"/api/jobs/{key}/download/site_truth.tsv")
                self.assertEqual(code, 200)
                truth = list(csv.DictReader(StringIO(data.decode()), delimiter="\t"))
                self.assertEqual([int(row["position"]) for row in truth], list(range(1, alignment.length + 1)))
                for index, row in enumerate(truth):
                    self.assertAlmostEqual(float(row["observed_gap_fraction"]),
                                           sum(seq[index] == "-" for seq in alignment.sequences) / size)
                if parameters.get("insertion_rate"):
                    self.assertGreater(alignment.length, ancestor)
                    self.assertEqual({row["mutation_rate_class"] for row in truth}, {"low", "medium", "high"})
                for method in ("rvne", "scaw-rvne"):
                    score_key, scored = self.job({
                        "filename": result["filename"], "fasta": result["fasta"], "mode": "aligned",
                        "method": method, "preset": "default", "maxGap": .3,
                    }, path="/api/jobs")
                    self.assertEqual(scored["state"], "done", scored)
                    self.assertEqual(len(scored["result"]["sites"]), alignment.length)
                    self.assertGreater(scored["result"]["scored"], 0)
                    self.assertEqual([int(site["position"]) for site in scored["result"]["sites"]],
                                     [int(row["position"]) for row in truth])
                    self.assertEqual(self.request(f"/api/jobs/{score_key}/download/scores.tsv")[0], 200)


if __name__ == "__main__":
    unittest.main()

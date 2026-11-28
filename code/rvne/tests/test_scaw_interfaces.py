"""Exercise the live server adapter and CLI, including exported effective settings."""
from contextlib import redirect_stderr
import importlib.util
from io import StringIO
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

from msa_conservation.cli import main

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code/frontend"))
spec = importlib.util.spec_from_file_location("frontend_server", ROOT / "code/frontend/server.py")
server = importlib.util.module_from_spec(spec)
spec.loader.exec_module(server)


class ScawInterfaceTests(unittest.TestCase):
    def test_live_presets_and_custom_download_match(self):
        for preset in ("default", "frozen", "custom"):
            body = {"example": "small-alignment", "mode": "aligned", "method": "scaw-rvne",
                    "aligner": "mafft", "preset": preset, "maxGap": .3}
            if preset == "custom":
                body["parameters"] = {"rvne": {"r": .8, "s1": .8, "s2": .8},
                                      "context": {"penalty_scale": 1.75}}
            options = server.load_options(body)
            alignment = server.read_fasta(StringIO(options["fasta"]))
            with TemporaryDirectory() as directory:
                result = server.run_scoring(alignment, options, Path(directory))
                effective = json.loads((Path(directory) / "settings.json").read_text())
                self.assertEqual(set(effective["context"]), {"penalty_scale"})
                self.assertEqual(effective["context_version"], "scaw-v2")
                self.assertEqual(effective["context_derived"]["solver"], "PELT")
                self.assertGreaterEqual(effective["context_derived"]["penalty"], 0)
                self.assertEqual(effective["context"]["penalty_scale"],
                                 1.75 if preset == "custom" else 1.)

    def test_stale_or_invalid_custom_context_fails(self):
        rvne = {"r": .8, "s1": .8, "s2": .8}
        for context in ({"beta": 6}, {"penalty_scale": True}, {"penalty_scale": -1},
                        {"penalty_scale": float("nan")}, {"penalty_scale": None}):
            with self.assertRaises(ValueError):
                server.custom_settings({"rvne": rvne, "context": context}, "scaw-rvne")

    def test_cli_override_and_v1_rejection(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)
            source = ROOT / "code/rvne/data/examples/small_alignment.fasta"
            settings, output = path / "settings.json", path / "sites.tsv"
            main(["score", "--input", str(source), "--method", "scaw-rvne",
                  "--penalty-scale", "2", "--settings-output", str(settings), "--output", str(output)])
            self.assertEqual(json.loads(settings.read_text())["context"], {"penalty_scale": 2.})
            obsolete = path / "obsolete.json"
            obsolete.write_text(json.dumps({"context": {"beta": 6, "smoothing": .25}}),
                                encoding="utf-8")
            error = StringIO()
            with redirect_stderr(error), self.assertRaises(SystemExit) as stopped:
                main(["score", "--input", str(source), "--method", "scaw-rvne",
                      "--settings", str(obsolete), "--output", str(output)])
            self.assertEqual(stopped.exception.code, 2)
            self.assertIn("v1 settings", error.getvalue())

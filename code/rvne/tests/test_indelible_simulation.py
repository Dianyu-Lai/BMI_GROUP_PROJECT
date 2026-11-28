from contextlib import redirect_stdout
import csv
from io import StringIO
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from benchmarks.run_supplemental_benchmark import load_families
from data_generation.cli import main
from data_generation.indelible import Options, find_executable
from data_generation.rvne_simulation import Region, scenario_rates, simulate, write_simulation


class RegionDesignTests(unittest.TestCase):
    def test_supplemental_regions_keep_their_original_layout(self):
        self.assertEqual(scenario_rates("conserved", 24, .4, 3), (.06,) * 24)
        self.assertEqual(scenario_rates("conserved-variable-blocks", 24, .25, 3),
                         (.06,) * 9 + (1.,) * 6 + (.06,) * 9)
        self.assertEqual(scenario_rates("single-motif", 24, .4, 3), (1.,) * 10 + (.06,) * 3 + (1.,) * 11)
        self.assertEqual(scenario_rates("three-rate-blocks", 25, .4, 3),
                         (.06,) * 4 + (.32,) * 3 + (1.,) * 10 + (.32,) * 3 + (.06,) * 5)

    def test_missing_executable_reports_how_to_configure_it(self):
        with patch("data_generation.indelible.shutil.which", return_value=None):
            with self.assertRaisesRegex(FileNotFoundError, "INDELIBLE_BIN"):
                find_executable("no-such-indelible-executable")


class IndelibleIntegrationTests(unittest.TestCase):
    def setUp(self):
        try:
            self.executable = find_executable()
        except FileNotFoundError:
            self.skipTest("INDELible 1.03 is required for executable integration tests")

    def test_redundant_family_with_short_tree_branches(self):
        result = simulate(size=50, length=240, scenario="redundant", seed=20107,
                          executable=self.executable)
        self.assertEqual((result.alignment.size, result.alignment.length), (50, 240))

    def test_cli_config_resolves_regions_and_applies_overrides(self):
        with TemporaryDirectory() as directory:
            folder = Path(directory)
            config = folder / "settings.json"
            config.write_text(json.dumps({"regions": "regions.json", "output_dir": "dataset",
                                          "size": 4, "seed": 42, "divergence": .9}), encoding="utf-8")
            (folder / "regions.json").write_text(json.dumps([
                {"length": 20, "rate": .06}, {"length": 20, "rate": .32},
                {"length": 20, "rate": 1., "insertion_rate": .08, "deletion_rate": .05},
            ]), encoding="utf-8")
            with redirect_stdout(StringIO()):
                self.assertEqual(main(["--config", str(config), "--seed", "43",
                                       "--indelible-bin", str(self.executable)]), 0)
            dataset = folder / "dataset"
            metadata = json.loads((dataset / "metadata.json").read_text(encoding="utf-8"))
            self.assertEqual((metadata["simulator"], metadata["seed"], metadata["ancestral_length"]),
                             ("INDELible", 43, 60))
            self.assertEqual(metadata["requested_divergence"], .9)
            self.assertEqual([part["rate"] for part in metadata["partitions"]], [.06, .32, 1.])
            self.assertEqual(metadata["partitions"][-1]["requested_options"]["insertion_rate"], .08)
            self.assertTrue((dataset / "alignment.fasta").is_file())
            with (dataset / "ground_truth.tsv").open() as handle:
                truth = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(truth), metadata["aligned_length"])
            self.assertEqual([int(row["position"]) for row in truth], list(range(1, len(truth) + 1)))

    def test_three_region_output_and_seed_reproducibility(self):
        regions = tuple(Region(20, rate) for rate in (.06, .32, 1.))
        settings = dict(size=4, length=60, regions=regions, executable=self.executable,
                        options=Options(insertion_rate=.08, deletion_rate=.05), boundary_width=2)
        with TemporaryDirectory() as directory:
            folder = Path(directory)
            first = simulate(**settings, seed=42, work_dir=folder / "first")
            second = simulate(**settings, seed=42, work_dir=folder / "second")
            altered = simulate(**settings, seed=43, work_dir=folder / "altered")
            self.assertEqual(first.alignment, second.alignment)
            self.assertNotEqual(first.alignment, altered.alignment)
            np.testing.assert_array_equal(first.rates, second.rates)
            np.testing.assert_allclose(first.truth, np.exp(-.7 * first.rates))
            self.assertEqual(first.alignment.size, 4)
            self.assertEqual(len(first.metadata["partitions"]), 3)
            self.assertIsNone(first.metadata["binary_metrics_reason"])
            write_simulation(first, folder / "family.fasta", folder / "truth.tsv", folder / "first")
            with (folder / "truth.tsv").open() as handle:
                truth = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(truth), first.alignment.length)
            self.assertEqual({row["mutation_rate_class"] for row in truth}, {"low", "medium", "high"})
            self.assertEqual(sum(int(row["boundary_mask"]) for row in truth), 8)
            for part in first.metadata["partitions"]:
                raw = folder / "first" / Path(part["control_file"]).parent
                for output in ("simulation_TRUE_1.fas", "simulation_ANCESTRAL_1.fas", "simulation_RATES.txt"):
                    self.assertTrue((raw / output).is_file())
            (folder / "manifest.tsv").write_text(
                "pattern\tseed\tsize\tlength\tfasta\ttruth\n"
                f"D\t42\t4\t{first.alignment.length}\tfamily.fasta\ttruth.tsv\n", encoding="ascii"
            )
            reloaded, = load_families(folder / "manifest.tsv")
            self.assertEqual(reloaded.truth.alignment, first.alignment)
            np.testing.assert_array_equal(reloaded.truth.conserved, first.conserved)

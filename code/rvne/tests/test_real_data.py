"""Check that the real-data CLI compares residues even when aligners shift gaps."""

from contextlib import redirect_stderr
import csv
from io import StringIO
import json
from pathlib import Path
from subprocess import CompletedProcess
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

from msa_conservation.cli import main


class RealDataTests(TestCase):
    def test_compare_by_reference_residue_and_keep_filtered_scores_blank(self):
        mafft = b">seq0000002\nAAD\n>seq0000001\nA-D\n>seq0000003\nA-D\n"
        clustalw = b">seq0000001\nAD-\n>seq0000002\nAAD\n>seq0000003\nA-D\n"

        def run(command, *, cwd, **kwargs):
            self.assertIn(">seq0000001\nAD\n", (cwd / "input.fasta").read_text(encoding="utf-8"))
            if Path(command[0]).name == "mafft":
                self.assertIn("--localpair", command)
                return CompletedProcess(command, 0, mafft, b"MAFFT v7.526")
            self.assertIn("-OUTPUT=FASTA", command)
            (cwd / "aligned.fasta").write_bytes(clustalw)
            return CompletedProcess(command, 0, b"CLUSTAL 2.1 Multiple Sequence Alignments", b"")

        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "raw.fasta"
            source.write_text(">ref|long-original-ID description\nAD\n>two\nAAD\n>three\nAD\n", encoding="utf-8")
            output = root / "results"
            with patch("msa_conservation.real_data.shutil.which", side_effect=lambda name: str(root / name)), \
                 patch("msa_conservation.real_data.subprocess.run", side_effect=run):
                self.assertEqual(main(["real-data", "--input", str(source), "--output-dir", str(output)]), 0)
            with (output / "reference_scores.tsv").open(encoding="utf-8") as handle:
                first, second = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(first["reference_position"], "1")
            self.assertEqual((second["reference_position"], second["reference_residue"]), ("2", "D"))
            self.assertEqual((second["mafft_column"], second["clustalw_column"]), ("3", "2"))
            self.assertTrue(second["mafft_score"])
            self.assertEqual(second["clustalw_score"], "")
            self.assertTrue((output / "mafft.fasta").read_text(encoding="utf-8").startswith(">ref|long-original-ID\n"))
            summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["aligners"]["mafft"]["version"], "7.526")
            self.assertEqual(summary["aligners"]["clustalw"]["version"], "2.1")
            self.assertEqual(summary["jointly_scored_reference_positions"], 1)
            self.assertIsNone(summary["reference_score_pearson"])

    def test_changed_residues_fail_before_publishing_outputs(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fasta"
            source.write_text(">ref\nAD\n>other\nAAD\n", encoding="utf-8")
            output = root / "results"
            error = StringIO()
            result = CompletedProcess(["mafft"], 0, b">seq0000001\nA-E\n>seq0000002\nAAD\n", b"")
            with patch("msa_conservation.real_data.shutil.which", return_value=str(root / "mafft")), \
                 patch("msa_conservation.real_data.subprocess.run", return_value=result), redirect_stderr(error):
                with self.assertRaises(SystemExit) as raised:
                    main(["real-data", "--aligner", "mafft", "--input", str(source), "--output-dir", str(output)])
            self.assertEqual(raised.exception.code, 2)
            self.assertIn("changed the ungapped sequence", error.getvalue())
            self.assertFalse(output.exists())

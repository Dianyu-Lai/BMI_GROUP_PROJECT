# BMI_GROUP_PROJECT

BMI3 Algorithms ICA Group 1 scores residue conservation in globular protein multiple sequence alignments. The main implementation reproduces the original relative von Neumann entropy (RVNE) score and evaluates segmentation-constrained adaptive windows (SCAW-RVNE) against published scoring methods on controlled synthetic alignments.

## Project structure

| Path | Contents |
| --- | --- |
| [`code/data_generation/`](code/data_generation/README.md) | Synthetic MSA generators and supplemental A–D datasets with ground truth. |
| [`code/benchmarks/`](code/benchmarks/README.md) | Parameter selection, held-out evaluation, ablations, timing and plots. |
| [`code/rvne/`](code/rvne/README.md) | Main implementation of RVNE, batched RVNE and SCAW-RVNE. |
| [`code/rvne/data/`](code/rvne/data/) | Example alignment, substitution-matrix inputs and the real-protein FASTA example. |
| [`code/rvne/results/reference_run/`](code/rvne/results/reference_run/) | Frozen parameters, benchmark tables and figures from the reference run. |
| [`code/rvne/tests/`](code/rvne/tests/) | Scoring, numerical equivalence, benchmark and CLI tests. |
| [`code/rvne/docs/`](code/rvne/docs/) | Method specifications, [benchmark results](code/rvne/docs/benchmark_report.md), [real-data protocol](code/rvne/docs/real_data_alignment.md) and [software documentation](code/rvne/docs/software_documentation.tex). |
| [`code/baselines/`](code/baselines/README.md) | Earlier entropy, BLOSUM62 and AAindex scoring implementations, CLI and tests. |

## Dependencies and setup

Both implementations require Python 3.10 or newer.

| Workflow | Additional dependencies |
| --- | --- |
| RVNE, SCAW-RVNE, synthetic benchmarks, tests and SVG plots | `numpy>=1.26,<3`, installed from [`code/rvne/requirements.txt`](code/rvne/requirements.txt). |
| Real-data alignment from unaligned FASTA | MAFFT 7 and ClustalW 2.1 for the default comparison. A single-aligner run requires the selected executable. See the [real-data setup](code/rvne/docs/real_data_alignment.md). |
| Earlier baseline experiments | Python standard library only. |

Run from the repository root in PowerShell to set up and verify the main implementation.

```powershell
cd code/rvne
python -m pip install -r requirements.txt
$env:PYTHONPATH = "..;src"
python -m unittest discover -s tests
```

In a POSIX shell, use `export PYTHONPATH=..:src` in place of the PowerShell assignment. Both scoring implementations expose the `msa_conservation` module name. See the [code guide](code/README.md) for the setup for each workflow.

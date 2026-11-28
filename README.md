# BMI_GROUP_PROJECT

BMI3 Algorithms ICA Group 1 scores residue conservation in globular protein multiple sequence alignments (MSAs). The project reproduces relative von Neumann entropy (RVNE), provides a batched implementation, and evaluates penalized segmentation with fixed boundary-clipped context (SCAW-RVNE v2) against published baselines on synthetic alignments with known conservation patterns. For unaligned protein FASTA, the real-data workflow can run MAFFT or ClustalW before scoring.

## Subprojects and key files

### RVNE, SCAW-RVNE, and real-data alignment

[Module guide](code/rvne/README.md)

| File or directory | Purpose |
| --- | --- |
| [`cli.py`](code/rvne/src/msa_conservation/cli.py) and [`rvne_cli.py`](code/rvne/src/msa_conservation/rvne_cli.py) | Dispatch the `score` and `real-data` commands. The scoring command writes per-site scores and optional segment and settings files. |
| [`msa.py`](code/rvne/src/msa_conservation/msa.py) | Parse and validate protein FASTA, calculate Henikoff sequence weights, and build column frequency and coverage profiles. |
| [`rvne.py`](code/rvne/src/msa_conservation/scoring/rvne.py) | Reproduce the original RVNE score and build the input profiles used by the scoring methods. |
| [`rvne_optimized.py`](code/rvne/src/msa_conservation/scoring/rvne_optimized.py) | Compute the same RVNE score with batched matrix operations and reused background calculations. |
| [`scaw.py`](code/rvne/src/msa_conservation/scoring/scaw.py) | Detect conservation regimes with PELT, clip a fixed window within each segment, and implement the smoothing and segmentation ablations. |
| [`substitution.py`](code/rvne/src/msa_conservation/scoring/substitution.py) | Build the BLOSUM62 substitution severity matrix used by the benchmark comparator. |
| [`real_data.py`](code/rvne/src/msa_conservation/real_data.py) | Call MAFFT and ClustalW, score their alignments, and map scores to the ungapped reference sequence. See the [real-data protocol](code/rvne/docs/real_data_alignment.md). |
| [`data/`](code/rvne/data/) and [`tests/`](code/rvne/tests/) | Example FASTA and reference probabilities, plus tests for scoring, numerical equivalence, benchmarking, and CLI integration. |
| [`docs/`](code/rvne/docs/) and [`results/`](code/rvne/results/) | Method specifications, complexity and benchmark reports, frozen parameters, and recorded results. |

### Published baselines

[Module guide](code/baselines/README.md)

| File or directory | Purpose |
| --- | --- |
| [`cli.py`](code/baselines/src/msa_conservation/cli.py) and [`msa.py`](code/baselines/src/msa_conservation/msa.py) | Read an aligned FASTA, optionally apply Henikoff weights, and export a per-site TSV. |
| [`baseline.py`](code/baselines/src/msa_conservation/scoring/baseline.py) | Assemble coverage and per-column scores. The default is Shannon entropy. `--reference-baselines` adds the four other scores. |
| [`entropy.py`](code/baselines/src/msa_conservation/scoring/entropy.py) | Implement normalized Shannon entropy, KL divergence, Jensen–Shannon divergence, and adjusted entropy. |
| [`substitution.py`](code/baselines/src/msa_conservation/scoring/substitution.py) and [`blosum62_joint.tsv`](code/baselines/data/reference/blosum62_joint.tsv) | Reconstruct integer BLOSUM62 scores and background frequencies from joint probabilities, and compute weighted sum-of-pairs scores. |
| [`tests/`](code/baselines/tests/) | Check MSA handling, score calculations, and baseline benchmark integration. |

### Synthetic data generation

INDELible generates aligned FASTA, per-column truth and metadata. The benchmark runners load those files or call the same simulator. The browser uses that simulator for generation and calls the scoring packages for analysis. [Module guide](code/data_generation/README.md)

| File | Purpose |
| --- | --- |
| [`indelible.py`](code/data_generation/indelible.py) | Configure and run INDELible 1.03 in isolated jobs, validate true alignments against unaligned sequences, and retain raw outputs and logs. |
| [`rvne_simulation.py`](code/data_generation/rvne_simulation.py) | Assemble region-specific evolutionary simulations on a shared rooted tree and export per-column truth, calibration measurements, and metadata. |
| [__main__.py](code/data_generation/__main__.py) and [cli.py](code/data_generation/cli.py) | Expose one CLI for individual families and benchmark suites with JSON settings and CLI overrides. |
| [suites.py](code/data_generation/suites.py) | Define patterns A–D, including calibrated gaps for D, and write family manifests. |

### Benchmarks and ablations

[Module guide](code/benchmarks/README.md)

| File | Purpose |
| --- | --- |
| [`run_benchmark.py`](code/benchmarks/run_benchmark.py) | Evaluate the standalone baseline scores on separate calibration and test alignments. [baseline_inputs.py](code/benchmarks/baseline_inputs.py) supplies INDELible rate inputs and matched-pair controls. |
| [`run_rvne_benchmark.py`](code/benchmarks/run_rvne_benchmark.py) | Tune RVNE and SCAW settings, evaluate held-out families, and run ablations, sensitivity sweeps, and timing comparisons. |
| [`run_supplemental_benchmark.py`](code/benchmarks/run_supplemental_benchmark.py) | Evaluate manifest-listed alignments and truth files using frozen parameters. |
| [`run_scaw_revision.py`](code/benchmarks/run_scaw_revision.py) | Run the SCAW v2 evaluation and validation sensitivity analysis with INDELible. |
| [`metrics.py`](code/benchmarks/metrics.py) and [`rvne_evaluation.py`](code/benchmarks/rvne_evaluation.py) | Calculate site and region metrics, select thresholds, match boundaries, and compute the main benchmark's reference comparator scores. |
| [`plot_rvne.py`](code/benchmarks/plot_rvne.py) | Generate SVG plots for the RVNE ablations and a short-motif example. |

### Browser demo

| File | Purpose |
| --- | --- |
| [`server.py`](code/frontend/server.py) | Serve the local interface and handle alignment, scoring, generation and downloads. |
| [`generation.py`](code/frontend/generation.py) | Validate simulator inputs and return aligned FASTA and site truth. |
| [`index.html`](code/frontend/index.html), [`app.js`](code/frontend/app.js) and [`results.js`](code/frontend/results.js) | Upload or select an MSA and inspect linked conservation and alignment views. |
| [`generate.js`](code/frontend/generate.js) and [`benchmarks.js`](code/frontend/benchmarks.js) | Run INDELible and display recorded benchmark results. |

## Dependencies

All Python workflows require **Python 3.10 or newer**. Install the project's Python dependency with [`code/rvne/requirements.txt`](code/rvne/requirements.txt).

| Workflow | Packages or tools |
| --- | --- |
| Standalone published baselines and matched-constraint control | Python standard library. |
| RVNE and SCAW-RVNE scoring, and supplemental evaluation | `numpy>=1.26,<3`. |
| Browser demo and its tests | NumPy and a browser with JavaScript module support. Simulator and aligner features use the tools below. |
| Synthetic evolutionary generation, baseline rate benchmark, RVNE reference benchmark, and synthetic motif plotting | NumPy and external **INDELible 1.03**. Use the project-local `.tools/indelible/` directory, `INDELIBLE_BIN`, or `PATH`. [Installation and dataset commands](code/data_generation/README.md). |
| Real-data alignment from unaligned FASTA | NumPy and external command-line tools **MAFFT 7** and **ClustalW 2.1** for the two-aligner comparison. A single-aligner run requires the selected executable. [Installation and executable paths](code/rvne/docs/real_data_alignment.md). |

## Quick start

Run from the repository root in PowerShell. Both scoring packages use the module name `msa_conservation`. Set `PYTHONPATH` for one scoring package at a time. Include `code` so Python can also find `benchmarks` and `data_generation`.

**RVNE and SCAW-RVNE**

```powershell
python -m pip install -r code/rvne/requirements.txt
$env:PYTHONPATH = "code;code/rvne/src"
python -m msa_conservation.cli score --method scaw-rvne --input code/rvne/data/examples/small_alignment.fasta
python -m unittest discover -s code/rvne/tests
```

**Published baselines**

```powershell
$env:PYTHONPATH = "code;code/baselines/src"
python -m msa_conservation.cli code/rvne/data/examples/small_alignment.fasta --reference-baselines
python -m unittest discover -s code/baselines/tests
```

For a POSIX shell, use `export PYTHONPATH=code:code/rvne/src` or `export PYTHONPATH=code:code/baselines/src` in place of the PowerShell assignment.

**Data generation, benchmarks and browser demo**

| To run | Command and guide |
| --- | --- |
| Check INDELible and set up a simulation | `python code/data_generation --check`. See the [generation guide](code/data_generation/README.md) for region and gap settings. |
| Tune and evaluate RVNE and SCAW | With the RVNE `PYTHONPATH`, run `python -m benchmarks.run_rvne_benchmark --output code/rvne/results/new_run`. See the [benchmark guide](code/benchmarks/README.md) and [recorded SCAW results](code/rvne/results/scaw_revision/README.md). |
| Start or test the browser demo | `python code/frontend/server.py` or `python -m unittest discover -s code/frontend/tests`. See the [browser guide](code/frontend/README.md). |

See the [SCAW method specification](code/rvne/docs/algorithm_design.md).

# BMI_GROUP_PROJECT

BMI3 Algorithms ICA Group 1 implements position-specific conservation scoring for aligned globular proteins. The project reproduces Johansson and Toh's relative von Neumann entropy (RVNE), tests an exact batched implementation, and evaluates a method that fits conservation regimes with penalized weighted SSE, then applies a fixed local window clipped to those regimes (SCAW v2).

## Installation

Run from `code/rvne` with Python 3.10 or newer. [`requirements.txt`](requirements.txt) declares `numpy>=1.26,<3` for scoring, benchmarking, SVG generation and tests. The benchmark was verified with Python 3.13.5 and NumPy 2.5.2. The `real-data` command also calls MAFFT 7 and/or ClustalW 2.1 through their command-line executables. See the [real-data setup](docs/real_data_alignment.md).

```powershell
python -m pip install -r requirements.txt
$env:PYTHONPATH = "..;src"
python -m unittest discover -s tests
```

In POSIX shells, use `export PYTHONPATH=..:src`.

## Score an alignment

```powershell
python -m msa_conservation.cli score --input data/examples/small_alignment.fasta --method rvne --output results/rvne.tsv
python -m msa_conservation.cli score --input data/examples/small_alignment.fasta --method scaw-rvne --output results/scaw.tsv --segments-output results/segments.tsv --settings-output results/scaw_settings.json
```

Use `rvne` for the original formula, `rvne-optimized` for the batched trace formula, and `uniform`, `gaussian`, `segmentation` or `scaw-rvne` for contextual methods. Omit `--output` to write the site table to stdout. The TSV gives raw base-20 divergence, nongap coverage, local estimate, final score, window size, support, variance and optional segment coordinates. Positions and segment endpoints are one-based and inclusive. Higher scores indicate stronger conservation. Sites with over 30% gaps or no canonical observations have blank scores.

The RVNE defaults are `--r 0.8 --s1 0.8 --s2 0.8 --max-gap 0.3`. SCAW v2 exposes `--penalty-scale` (default 1), which scales the penalty for each changepoint. It uses a fixed radius of 3 and center weight of 1/2, with windows clipped at segment boundaries. See the [penalty formula, assumptions and citations](docs/algorithm_design.md).

Use `--settings results/scaw_revision/parameters.json` for the v2 study configuration. Downloads record the method version, one context parameter, effective penalty and fixed design constants. SCAW settings files use the v2 schema with `penalty_scale`. Optional `--segment-thresholds LOWER UPPER` labels exported segment means, and needs independent calibration; otherwise classes are unclassified. These reporting cutoffs are separate from the one context parameter.

The [current benchmark](results/scaw_revision/README.md) evaluates SCAW and two ablations over 104 families, including 12 fresh INDELible simulations.

## Align and score real proteins

Install MAFFT 7 and ClustalW 2.1, then run from `code/rvne` with the same `PYTHONPATH` as above. The example aligns five UniProt globins with both tools and compares SCAW-RVNE scores in human myoglobin coordinates.

```powershell
python -m msa_conservation.cli real-data --input data/real/globins.fasta --reference-id P02144 --output-dir results/globins
```

See the [real-data protocol](docs/real_data_alignment.md) for the comparison design, Windows setup, input rules and output files.

## Evaluate RVNE on synthetic data

Use the [data generators](../data_generation/README.md) to create aligned FASTA files and truth labels, and the [benchmark commands](../benchmarks/README.md) to tune and evaluate RVNE or reproduce the ablations.

## Repository guide

| Path | Contents |
| --- | --- |
| `src/msa_conservation/` | FASTA parsing, sequence weighting and CLI entry points |
| `src/msa_conservation/scoring/` | Original RVNE, batched RVNE, contextual scoring and SCAW |
| `src/msa_conservation/real_data.py` | MAFFT and ClustalW integration, scoring and reference residue mapping |
| [`../benchmarks/`](../benchmarks/README.md) | Parameter selection, evaluation, ablations, timing and SVG plots |
| `data/examples/` | Example aligned FASTA for scoring |
| `data/real/` | Unaligned UniProt globin sequences for the real-data workflow |
| `data/reference/` | BLOSUM62 pair probabilities for the benchmark comparator |
| [`../data_generation/datasets/`](../data_generation/README.md) | Generated alignments and truth labels |
| `tests` | Mathematical checks, optimization equivalence and CLI integration |
| `results/scaw_revision` | Current SCAW benchmark, simulation inputs and v2 settings |
| `docs` | Methods, complexity, benchmark protocol, results and LaTeX documentation |

- [Original RVNE specification](docs/rvne_original_specification.md)
- [Algorithm design](docs/algorithm_design.md)
- [Computational complexity](docs/complexity_analysis.md)
- [Benchmark protocol](docs/benchmark_protocol.md) and [results](docs/benchmark_report.md)
- [ICA software documentation](docs/software_documentation.tex)

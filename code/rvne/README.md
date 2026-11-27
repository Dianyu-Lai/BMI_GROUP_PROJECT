# BMI_GROUP_PROJECT

BMI3 Algorithms ICA Group 1 implements position-specific conservation scoring for aligned globular proteins. The project reproduces Johansson and Toh's relative von Neumann entropy (RVNE), tests an exact batched implementation, and evaluates a method that fits conservation regimes before selecting a local smoothing window.

## Installation

Run from `code/rvne` in the BMI3 repository with Python 3.10 or newer. The benchmark was verified with Python 3.13.5 and NumPy 2.5.2.

```powershell
python -m pip install -r requirements.txt
$env:PYTHONPATH = "src"
python -m unittest discover -s tests
```

In POSIX shells, use `export PYTHONPATH=src`.

## Score an alignment

```powershell
python -m msa_conservation.cli score --input data/examples/small_alignment.fasta --method rvne --output results/rvne.tsv
python -m msa_conservation.cli score --input data/examples/small_alignment.fasta --method scaw-rvne --output results/scaw.tsv --segments-output results/segments.tsv --settings-output results/scaw_settings.json
```

Use `rvne` for the original formula, `rvne-optimized` for the batched trace formula, and `uniform`, `gaussian`, `adaptive`, `segmentation` or `scaw-rvne` for contextual methods. Omit `--output` to write the site table to stdout. The TSV gives raw base-20 divergence, nongap coverage, local estimate, final score, window size, support, variance and optional segment coordinates. Positions and segment endpoints are one-based and inclusive. Higher scores indicate stronger conservation. Sites with over 30% gaps or no canonical observations have blank scores. Example site 3 has raw score 1.3219 and default SCAW score 1.2998 from a three-column window.

The scoring defaults are `--r 0.8 --s1 0.8 --s2 0.8 --max-gap 0.3`. The context defaults are `--beta 6 --smoothing 0.25 --radii 0 1 2 4 8 --min-support 1.5 --min-segment 3` and an unbounded maximum length. Each flag is optional. `--r` controls background interpolation, `--s1` and `--s2` control observed and background similarity. `--beta` sets the segmentation penalty in global-variance units. `--smoothing` sets the weight on reciprocal support in the window objective. `--max-segment` can constrain segment length and takes zero for an unbounded partition. Export raw-mean segment classes using two increasing, independently calibrated `--segment-thresholds`. The default segment class is unclassified. Use `--settings results/reference_run/parameters.json` to select the benchmark's retuned configuration.

## Simulate and benchmark

```powershell
python -m benchmarks.rvne_simulation --scenario motifs --seed 42 --output data/synthetic/motifs.fasta --truth data/synthetic/motifs_truth.tsv
$env:OPENBLAS_NUM_THREADS = "1"
python -m benchmarks.run_rvne_benchmark --mode all --sweep --repeats 4 --output results/rvne
python -m benchmarks.plot_rvne --results results/rvne
```

`--mode tune` writes the tuning grid and frozen settings without evaluating the test sets. `--mode benchmark --parameters results/reference_run/parameters.json` evaluates the held-out data using the saved settings. Wall times vary across runs. Keep `results/reference_run` as the recorded run and direct new runs to another directory.

### Generate supplemental test datasets

```powershell
python -m benchmarks.generate_pattern_data --output data/synthetic/patterns_abc
```

Choose a new output directory. The generator writes 36 aligned FASTA files, matching site-truth TSVs, and a manifest with the settings and seeds. Each alignment has 50 sequences and 240 columns. Each setting has four independent replicates.

| Pattern | Site-rate profile | Varied setting |
| --- | --- | --- |
| A | Uniformly conserved at rate 0.06 | Divergence 0.35, 0.7, 1.4 |
| B | Conserved flanks at 0.06 and a central variable block at 1.0 | Variable block fraction 25%, 40%, 60% |
| C | A central conserved motif at 0.06 in a variable background at 1.0 | Motif length 3, 5, 9 residues |

The truth TSV records the one-based position, simulated mutation rate, expected conservation probability, conserved label and a flag for a boundary immediately before the position. The supplemental seeds 50100–50135 are independent of the reference benchmark's seeds.

Generate the optional three-rate Pattern D in a new directory.

```powershell
python -m benchmarks.generate_pattern_data --suite d --output data/synthetic/patterns_d
```

Pattern D uses five consecutive blocks with low, medium, high, medium and low mutation rates (0.06, 0.32 and 1.0). The central high-rate block takes up 25%, 40% or 60% of the alignment, matching Pattern B's settings. The medium-rate blocks take up half of each flank, rounded down to a whole number of columns. At divergence 0.7, the three per-descendant mutation probabilities are approximately 4.1%, 20.1% and 50.3%, using `1 - exp(-divergence * rate)`. The suite contains 12 families (four seeds per setting, 50136–50147), with 50 sequences and 240 columns each.

Use Spearman correlation against the `conservation` column to check whether scores recover the three-tier ordering. The binary `conserved` column marks only the low-rate sites as conserved. The truth file records all four rate transitions for boundary evaluation.

## Repository guide

| Path | Contents |
| --- | --- |
| `src/msa_conservation` | FASTA parsing, sequence weighting, RVNE, exact optimization, contextual scoring and CLI |
| `benchmarks` | Synthetic families, metrics, parameter selection, ablation and timing runners |
| `data/examples` | Runnable example alignment |
| `tests` | Mathematical checks, optimization equivalence and CLI integration |
| `results/reference_run` | Frozen settings, held-out tables, timing and figures |
| `docs` | Methods, complexity, benchmark protocol, results and LaTeX documentation |

- [Original RVNE specification](docs/rvne_original_specification.md)
- [Algorithm design](docs/algorithm_design.md)
- [Computational complexity](docs/complexity_analysis.md)
- [Benchmark protocol](docs/benchmark_protocol.md) and [results](docs/benchmark_report.md)
- [ICA software documentation](docs/software_documentation.tex)

OpenAI Codex assisted the implementation, testing, benchmarking and documentation. The previous exploratory AAindex/BLOSUM pipeline and its commands are preserved in [the earlier project notes](docs/earlier_pipeline.md).

# Conservation benchmarks

Run the RVNE benchmarks from `code/` in PowerShell. Install [the RVNE requirements](../rvne/requirements.txt) before running these commands.

```powershell
$env:PYTHONPATH = "rvne/src"
$env:OPENBLAS_NUM_THREADS = "1"
python -m benchmarks.run_rvne_benchmark --mode all --sweep --repeats 4 --output rvne/results/rvne
python -m benchmarks.plot_rvne --results rvne/results/rvne
```

The generator uses INDELible 1.03. The substitution-model sweep compares WAG, JTT and LG.

In a POSIX shell, use `export PYTHONPATH=rvne/src` and `export OPENBLAS_NUM_THREADS=1`.

`--mode tune` writes the tuning grid and frozen settings. Use `--mode benchmark --parameters rvne/results/scaw_revision/parameters.json` to evaluate the held-out families with the recorded settings. Write new results to a separate directory to preserve the recorded run.

Generate the supplemental A–D datasets using the [data generation instructions](../data_generation/README.md), then pass their manifests to the evaluator.

```powershell
python -m benchmarks.run_supplemental_benchmark --manifest data_generation/datasets/indelible_abc/manifest.tsv --manifest data_generation/datasets/indelible_d/manifest.tsv --parameters rvne/results/rvne/parameters.json --output rvne/results/indelible_supplemental
```

| Module | Purpose |
| --- | --- |
| `run_rvne_benchmark.py` | Parameter selection, held-out evaluation, ablations, sensitivity sweeps and timing. |
| `run_supplemental_benchmark.py` | Frozen-parameter evaluation of supplied FASTA and site-truth pairs. |
| `run_benchmark.py` | [Published entropy and BLOSUM62 baseline comparisons](../baselines/README.md). Uses the scoring package under `baselines/src`. |
| `metrics.py`, `rvne_evaluation.py` | Evaluation metrics and reference comparator scores. |
| `plot_rvne.py` | SVGs for the RVNE ablations and a short-motif example. |

## Run the SCAW evaluation with INDELible

From the repository root:

```powershell
$env:PYTHONPATH = "code/rvne/src;code"
$env:OPENBLAS_NUM_THREADS = "1"
python -m benchmarks.run_scaw_revision --output code/rvne/results/scaw_revision_repeat --fresh-indelible
```

Use a new output directory. The run evaluates 32 INDELible core families with thresholds calibrated on 24 separate INDELible validation families. The command also includes the 12 dedicated INDELible cases. Add a repeatable `--manifest` for any generated supplemental suite. The penalty multiplier, radius and center weight remain fixed at 1, 3 and 1/2. For the original cohort and threshold provenance, see the [recorded study](../rvne/results/scaw_revision/README.md).

The general tuning runner searches c in {1, 0.5, 2} after RVNE selection and writes the context version. For supplemental scoring, use v2 settings and thresholds from that runner. The SCAW evaluation above uses its own threshold names.

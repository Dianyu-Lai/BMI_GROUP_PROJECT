# Conservation benchmarks

Run the RVNE benchmarks from `code/` in PowerShell. Install [the RVNE requirements](../rvne/requirements.txt) before running these commands.

```powershell
$env:PYTHONPATH = "rvne/src"
$env:OPENBLAS_NUM_THREADS = "1"
python -m benchmarks.run_rvne_benchmark --mode all --sweep --repeats 4 --output rvne/results/rvne
python -m benchmarks.plot_rvne --results rvne/results/rvne
```

The default generator uses INDELible 1.03. To reproduce the historical reference benchmark and the motif example, pass `--legacy-reference` to both commands. The current one-factor substitution-model sweep compares WAG with JTT and LG. The historical sweep uses its original conservative and radical substitution rules.

In a POSIX shell, use `export PYTHONPATH=rvne/src` and `export OPENBLAS_NUM_THREADS=1`.

`--mode tune` writes the tuning grid and frozen settings. Use `--mode benchmark --parameters rvne/results/reference_run/parameters.json` to evaluate the held-out families with the recorded settings. Write new results to a separate directory to preserve the recorded run.

Generate the supplemental A–D datasets using the [data generation instructions](../data_generation/README.md), then pass both manifests to the evaluator. Historical manifests from the ZIP archives are also accepted.

```powershell
python -m benchmarks.run_supplemental_benchmark --manifest data_generation/datasets/indelible_abc/manifest.tsv --manifest data_generation/datasets/indelible_d/manifest.tsv --parameters rvne/results/reference_run/parameters.json --output rvne/results/indelible_supplemental
```

| Module | Purpose |
| --- | --- |
| `run_rvne_benchmark.py` | Parameter selection, held-out evaluation, ablations, sensitivity sweeps and timing. |
| `run_supplemental_benchmark.py` | Frozen-parameter evaluation of supplied FASTA and site-truth pairs. |
| `run_benchmark.py` | [Published entropy and BLOSUM62 baseline comparisons](../baselines/README.md). Uses the scoring package under `baselines/src`. |
| `metrics.py`, `rvne_evaluation.py` | Evaluation metrics and reference comparator scores. |
| `plot_rvne.py` | SVGs for the RVNE ablations and a short-motif example. |

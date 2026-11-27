# BMI3 conservation scoring code

| Directory | Contents |
| --- | --- |
| [`data_generation/`](data_generation/README.md) | Synthetic MSA generators, A–D datasets, truth labels and manifests. |
| [`benchmarks/`](benchmarks/README.md) | Scoring comparisons, parameter tuning, held-out evaluation, ablations, metrics and plots. |
| [`rvne/`](rvne/README.md) | Original RVNE, batched RVNE, SCAW, and [real-protein alignment with MAFFT and ClustalW](rvne/docs/real_data_alignment.md). |
| [`baselines/`](baselines/README.md) | Earlier Shannon, KL, Jensen–Shannon, BLOSUM62 and adjusted-entropy scoring implementations. |
| [`output/`](output/benchmark_analysis.md) | Analysis of the recorded benchmark results, including tables and figures. |

Both scoring packages use the `msa_conservation` name. Follow the workflow README to select one scoring package on `PYTHONPATH`.

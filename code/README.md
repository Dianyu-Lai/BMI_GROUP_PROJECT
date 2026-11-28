# BMI3 conservation scoring code

| Directory | Contents |
| --- | --- |
| [`data_generation/`](data_generation/README.md) | One INDELible 1.03 command for individual protein MSAs and A–D suites with site truth, gap calibration and raw run metadata. |
| [`benchmarks/`](benchmarks/README.md) | Scoring comparisons, parameter tuning, held-out evaluation, ablations, metrics and plots. |
| [`rvne/`](rvne/README.md) | Original RVNE, batched RVNE, SCAW, and [real-protein alignment with MAFFT and ClustalW](rvne/docs/real_data_alignment.md). |
| [`baselines/`](baselines/README.md) | Shannon, KL, Jensen–Shannon, BLOSUM62 sum-of-pairs and adjusted-entropy scoring implementations. |
| [SCAW benchmark results](rvne/results/scaw_revision/README.md) | Current method and ablation results, parameters, and input fingerprints. |
| [`frontend/`](frontend/README.md) | Local FASTA upload, alignment, scoring, INDELible generation and benchmark plots. |

Both scoring packages use the `msa_conservation` name. Follow the workflow README to select one scoring package on `PYTHONPATH`.

Python 3.10+ and [`rvne/requirements.txt`](rvne/requirements.txt) cover the Python side of evolutionary data generation. Install **INDELible 1.03** separately and set `INDELIBLE_BIN` to its executable, or add it to `PATH`. [Setup and CLI examples](data_generation/README.md).

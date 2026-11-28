# SCAW v2 benchmark

Run date: 2026-11-28. [Method specification](../../docs/algorithm_design.md), [parameters](parameters.json) and [input fingerprints](inputs.tsv).

The demo uses the 12 INDELible families and reports metrics that do not require the historical decision thresholds. The 32 core families, 60 supplemental families and 24 threshold-calibration families all used the earlier independent-descendant simulator. Their outputs below are archived for reproduction.

## Design fixed before evaluation

The three scoring variants receive the same raw RVNE profile, using `(r,s1,s2)=(0,0.1,0)`, the same mask and coverage. SCAW uses a fixed penalty multiplier of 1, radius 3 and center weight 1/2. The validation sensitivity tables also report multipliers 0.5 and 2.

Each binary threshold was freshly calibrated for its own method on the 24 historical validation families. No threshold was fit to the test or fresh INDELible labels. The first two cohorts reuse previously examined data and are **not** newly independent test evidence. The 12 fresh INDELible 1.03 families use seeds 81101–81106 and 81201–81206, 20 sequences, 120 ancestral columns, WAG on the simulator tree, and homogeneous, variable-block, 3/5-column motif, three-rate and indel cases. Indel alignments can have different final lengths; exact input hashes and dimensions are recorded.

`summary.tsv` contains every method and scenario. Correlations and curve metrics exclude undefined single-class/constant-truth cases and report their own n. Boundary F1 uses one-to-one matching within 2 alignment columns; exact-boundary F1 is also recorded. Gradual profiles have no discrete boundary target and are excluded; no-change families are included (F1=1 only if no boundary is predicted). Region IoU thresholds entire segment means times their mean usable coverage, matching the historical diagnostic; it is not a one-to-one segment overlap score. Over/under-segmentation counts are numbers of extra/missing boundaries, not rates. Raw per-family data supports alternate summaries.

## Results by cohort

### 32 reused historical core families

| Method | Spearman | AUPRC | Boundary F1 ±2 | Region IoU | Specificity |
|---|---|---|---|---|---|
| v2 | 0.8409 (n=24) | 0.9970 (n=24) | 0.9898 (n=28) | 0.9279 (n=28) | 0.9474 (n=28) |
| v2_segment_mean | 0.8685 (n=24) | 0.9949 (n=24) | 0.9898 (n=28) | 0.9279 (n=28) | 0.9478 (n=28) |
| fixed_uniform_7 | 0.8428 (n=24) | 0.9982 (n=24) | N/A | N/A | 0.9501 (n=28) |

### 60 archived supplemental families

| Method | Spearman | AUPRC | Boundary F1 ±2 | Region IoU | Specificity |
|---|---|---|---|---|---|
| v2 | 0.7285 (n=48) | 0.9999 (n=48) | 0.9378 (n=60) | 0.9875 (n=60) | 0.9730 (n=48) |
| v2_segment_mean | 0.7705 (n=48) | 0.9992 (n=48) | 0.9378 (n=60) | 0.9982 (n=60) | 0.9993 (n=48) |
| fixed_uniform_7 | 0.7284 (n=48) | 0.9656 (n=48) | N/A | N/A | 0.9787 (n=48) |

### 12 fresh INDELible families

| Method | Spearman | AUPRC | Boundary F1 ±2 | Region IoU | Specificity |
|---|---|---|---|---|---|
| v2 | 0.6463 (n=10) | 0.8986 (n=10) | 0.4921 (n=12) | 0.5731 (n=12) | 0.7054 (n=10) |
| v2_segment_mean | 0.7546 (n=8) | 0.6099 (n=10) | 0.4921 (n=12) | 0.5731 (n=12) | 0.8144 (n=10) |
| fixed_uniform_7 | 0.6553 (n=10) | 0.8932 (n=10) | N/A | N/A | 0.7518 (n=10) |

On the 12 fresh INDELible families, SCAW had mean AUPRC 0.8986 and boundary F1 0.4921. The uniform-window ablation reached AUPRC 0.8932. These results quantify performance on simulated conservation labels.

## Context-only timing

Three measured repetitions after a warm-up on a noisy four-step profile with one OpenBLAS thread. Timings cover only the SCAW context stage.

| Columns | Median milliseconds |
|---|---:|
| 240 | 6.41 |
| 1000 | 28.35 |
| 5000 | 190.98 |

## Reproduction and outputs

The tables in this directory record the original evaluation. To rerun the study design on INDELible inputs, use the [current SCAW runner](../../../benchmarks/README.md#run-the-scaw-evaluation-with-indelible). Its 24 validation families supply newly calibrated thresholds, and its 32 core families and optional supplemental suites use INDELible sequences.

The local fresh input directory retains simulator controls, true alignments, rates, logs, derived FASTA and site truth. The `fresh_inputs.zip` archive contains the same 12 inputs; Git ignores the expanded directory. `inputs.tsv` fingerprints the aligned sequence bodies and rates. `parameters.json` records hashes of the method and runner used in the original study. `families.tsv`, `summary.tsv`, `validation_sensitivity.tsv`, `validation_summary.tsv` and `runtime.tsv` retain the original measurements. The frontend uses this directory for the study preset, examples and benchmark panel.

`parameters.json` contains method-named thresholds for this run. It is valid for scoring, but the general supplemental runner needs thresholds produced by a new generic tuning run; these two threshold namespaces must not be mixed.

### Matched column-only baselines

`baseline_summary.tsv` compares Original RVNE at the study parameters, Shannon, KL, Jensen–Shannon, BLOSUM62 sum of pairs, adjusted entropy, BLOSUM62 substitution cost, ClustalW consensus and AL2CO variance on the same 12 fresh families. The [published baselines](../../../baselines/README.md) and [RVNE protocol](../../docs/benchmark_protocol.md#comparators-and-metrics) define the scores. The two BLOSUM62 methods use different calculations and are reported separately.

All methods use SCAW's column mask and nongap coverage multiplier. The shared RVNE benchmark supplies Shannon, KL, Jensen–Shannon, BLOSUM62 substitution cost, ClustalW consensus and AL2CO variance. The standalone package supplies the distinct-pair BLOSUM62 mean and adjusted entropy using Henikoff weights. Before coverage multiplication, the script maps the BLOSUM62 mean to [0, 1] using the minimum and maximum integer matrix entries, and maps adjusted entropy to conservation as 1 / (1 + adjusted entropy). The script checks input fingerprints and the scored-column mask. It reports threshold-independent site metrics. The frontend displays AUPRC and Spearman before the window and segmentation ablations.

```powershell
$env:PYTHONPATH = "code/rvne/src;code"
python -m benchmarks.compare_scaw_baselines --results code/rvne/results/scaw_revision
```

## Verification and limits

Mathematical tests compare PELT to exhaustive weighted partitions over 192 randomized short profiles, test zero penalty/constant ties, missing profiles, adjacent-only noise estimation, fixed blending, short motif boundaries and affine transformations. Interface tests exercise CLI flags, settings versioning, stale-configuration rejection and the frontend adapter's default/study/custom scoring paths and exported settings. The general RVNE tuning runner was also executed through parameter selection and threshold calibration. See `verification.json` for the final checks.

Not yet established: empirical inverse-coverage variance, independent/Gaussian RVNE residuals, noise-estimator robustness to phylogeny and serial correlation, window constants optimal for RVNE, comprehensive sensitivity under alignment error, significance/false-discovery calibration, real protein functional-label performance or broad phylogenetic generalization. No full replacement of the historical large sensitivity sweep was run. Only local changes were made; nothing was pushed.

The method specification describes the PELT objective, fixed-window heuristic and first-difference estimator adaptations. [Killick et al.](https://doi.org/10.1080/01621459.2012.737745) · [Capra & Singh](https://doi.org/10.1093/bioinformatics/btm270) · [Dai et al.](https://doi.org/10.1214/17-STS613)

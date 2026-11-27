# Supplemental A/B/C/D benchmark

This run evaluates the supplied aligned FASTA and site-truth TSV pairs for 48 synthetic families. Each pattern has 12 families with 50 sequences and 240 columns per alignment. The runner evaluates E0–E8 and six column-level baselines using the original frozen parameters and method-specific validation thresholds. All columns received scores. Each result table aggregates equally weighted families and records the number contributing to each metric.

## Main results

| Pattern | Main observation |
| --- | --- |
| A, uniformly conserved | E8 made zero variable-site calls at the frozen threshold, while its segmentation introduced boundaries in 7 of 12 families, averaging 1.33 false boundaries per family. E5 and E6 averaged 0.67. AUROC and Spearman are undefined because the truth is constant within each family. |
| B, central variable block | E5, E6 and E8 achieved 1.000 variable-site IoU and 1.000 boundary F1. E7 and Shannon already had 1.000 variable-site recall, so this pattern mainly verifies that the block and its endpoints are recovered. |
| C, one short conserved motif | E8 recovered every motif site at widths of 3, 5 and 9 residues, with mean boundary F1 of 1.000. E7 and Shannon also recovered every motif site. E2 uniform and E3 Gaussian smoothing recovered none of the 3- or 5-residue motif sites under their frozen thresholds. |
| D, three mutation-rate levels | E8 detected 98.76% of medium- and high-rate sites pooled within each family, up from 80.55% with E7 and Shannon. The E8 recall improvement over E7 occurred in all 12 families. The gain came from medium-rate sites. Both E7 and E8 recovered 100% of high-rate sites. |

## D ablation and baselines

| Method | Spearman | Variable recall | Medium called variable | Low vs medium AUROC | Medium vs high AUROC | Boundary F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| E0 original RVNE | 0.8991 | 0.7061 | 0.3029 | 0.9554 | 0.9959 | — |
| E4 adaptive windows | 0.9158 | 0.6573 | 0.1836 | 0.9968 | 0.9997 | — |
| E5 segmentation only | **0.9743** | 0.7763 | 0.4572 | 0.9968 | 1.0000 | 1.0000 |
| E6 default SCAW-RVNE | 0.9154 | 0.6486 | 0.1590 | 0.9958 | 0.9997 | 1.0000 |
| E7 retuned RVNE | 0.9165 | 0.8055 | 0.5387 | 0.9963 | 0.9996 | — |
| E8 retuned SCAW-RVNE | 0.9175 | **0.9876** | **0.9530** | 0.9989 | 1.0000 | 1.0000 |
| Shannon | 0.9165 | 0.8055 | 0.5387 | 0.9963 | 0.9996 | — |
| AL2CO variance | 0.9161 | 0.8271 | 0.5873 | 0.9960 | 0.9995 | — |

D uses low, medium and high mutation rates of 0.06, 0.32 and 1.00. The frozen binary decision treats the medium and high rates as variable, and a higher score means greater conservation. The adjacent-rate AUROCs measure ranking with the lower rate as the positive class. The segmentation boundary metric counts a one-to-one match within two aligned columns and does not require the neighboring regimes to receive different binary labels. The separate site IoU uses thresholded column scores, while `region_iou` uses thresholded segment means and conserved-region ground truth.

E7 and E8 use the same retuned RVNE profiles. Their comparison measures the combined effect of contextual processing and the validation threshold for each output scale. E5 had the strongest three-level rank correlation among the methods evaluated. E5, E6 and E8 matched D's true boundaries within the two-column tolerance. E6's lower medium-rate recall shows that the default adaptive windows and segmentation do not by themselves account for E8's result. These alignments have no simulated gaps and share the simulator used in the original benchmark. The results support the specified synthetic profiles.

## Reproduce the run

From `code/`, install `rvne/requirements.txt` and use PowerShell. Use the archived A/B/C and D manifests under `data_generation/datasets/` and choose a new output directory.

```powershell
$env:PYTHONPATH = 'rvne/src'
$env:OPENBLAS_NUM_THREADS = '1'
python -m benchmarks.run_supplemental_benchmark --manifest data_generation/datasets/patterns_abc/manifest.tsv --manifest data_generation/datasets/patterns_d/manifest.tsv --parameters rvne/results/reference_run/parameters.json --output rvne/results/supplemental_repro
```

The three TSV outputs were identical when rerun from regenerated inputs. The original 480-row held-out table also reproduced exactly.

- [Per-family metrics](families.tsv)
- [Summary by pattern and method](summary.tsv)
- [Summary by divergence, block fraction or motif width](settings_summary.tsv)
- [Run metadata](run_metadata.json)

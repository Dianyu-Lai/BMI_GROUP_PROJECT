# RVNE benchmark protocol

See the [recorded SCAW evaluation](../results/scaw_revision/README.md) for the fixed-parameter run and its cohorts. The general runner also supports tuning and held-out evaluation.

## Ground truth and families

The general benchmark uses INDELible 1.03 with a shared rooted tree and empirical amino-acid substitution models. Designed site mutation rates supply the ground truth.

The core scenarios are uniformly conserved, uniformly variable, mixed blocks, short motifs, gradual transitions, high divergence, indels and redundancy. Blocks contain rates 0.06, 0.32 and 1.00. Motifs use 0.06 against 1.00. The gradual scenario approximates a cosine profile from 0.06 to 1.00 with up to 12 constant-rate partitions. The positive class is t_i ≤ 0.18. Core defaults are 50 sequences, 240 ancestral residues, root-to-tip divergence 0.7, variable fraction 40% and five-column motifs. High divergence multiplies divergence by 2.5. The indel scenario calibrates to 20% target gap occupancy. The redundancy scenario uses a near-tip clade with a nominal fraction of 70%. One-factor sweeps vary depth, length, divergence, motif width, variable fraction, gap occupancy, redundancy and the substitution model (WAG, JTT or LG). The runner records the exact levels.

## Split and selection

Family seeds are 100 × seed group plus settings index. Training uses groups 101–103, validation 201–203, core test 301–304 and sweeps 401–403. Core sets have 24 training, 24 validation and 32 test families. The sweep has 63 additional families.

Training ranks all 1,331 paper settings by mean family Spearman. Keep the 12 best, the paper default and the clipped ±0.05 Cartesian neighborhood of the top training setting. Validation selects among these on the same metric. At ties to 12 decimals, prefer settings nearer the paper reference in L1 distance, then lexicographic order. Constant-truth families do not contribute correlations. After freezing RVNE, select the SCAW penalty multiplier from 1, 0.5 and 2 by validation boundary F1. Ties prefer the first candidate. The study results instead fix the multiplier at 1 before test evaluation.

## Comparators and metrics

| ID | Method |
| --- | --- |
| E0 | Original RVNE, direct |
| E1 | Original RVNE, cached trace implementation |
| E2 | Original RVNE plus uniform smoothing |
| E3 | Original RVNE plus Gaussian smoothing |
| E5 | Original RVNE plus segment mean |
| E6 | Original RVNE plus fixed windows clipped to fitted boundaries |
| E7 | Retuned RVNE |
| E8 | Retuned RVNE plus tuned fixed windows clipped to fitted boundaries |

The uniform and Gaussian windows use radius 3. SCAW uses the same radius and clips the window to the fitted segment. All comparators use the RVNE site mask, nongap coverage multiplier and high-conservation orientation.

E4 is the historical adaptive-window ablation. Its implementation and results are available at commit `a4f94ee` under `code/rvne/src/msa_conservation/scoring/scaw.py` and `code/rvne/results/reference_run/`. The current runner uses the v2 method list above.

| Comparator | Per-column score before coverage weighting |
| --- | --- |
| Shannon | One minus base-20 entropy of the Henikoff frequencies |
| KL | Base-2 KL divergence of the Henikoff frequencies from the BLOSUM62 `.qij` marginals |
| JSD | Equal-weight Jensen–Shannon divergence in bits using the same distributions |
| BLOSUM62 | One minus the expected pair substitution cost divided by its matrix maximum, using the Henikoff frequencies |
| ClustalW consensus | [Published consensus categories](https://www.genome.jp/tools/clustalw/clustalw_readme.html) applied to the supplied alignment. Exact, strong-group, weak-group and blank columns use scores 3, 2, 1 and 0. Gaps and ambiguous residues give the blank category. |
| AL2CO variance | [Pei and Grishin's unweighted variance measure (2_1)](https://doi.org/10.1093/bioinformatics/17.8.700). Sum the squared differences between the column's canonical-residue frequencies and the canonical-residue frequencies pooled across the input alignment. |

The ClustalW comparison measures its consensus annotation rule on the same MSA as the other methods. The AL2CO comparison uses the raw single-column index before optional family-level normalization. These comparators implement the documented scoring rules in Python.

AUROC uses average ranks. AUPRC uses noninterpolated average precision with ties grouped. Pearson and Spearman are per family and undefined for a constant target or score. Spearman and tie-sensitive classification use scores rounded to 12 decimals. Sensitivity and specificity use each method's threshold, selected to maximize balanced accuracy on pooled validation sites. The general runner recalibrates thresholds on its validation split. Single-class families contribute whichever metric is defined. Report the scorable fraction and restrict position metrics to shared scorable columns. Report means, sample standard deviations and metric-specific counts across families.

Score abrupt boundaries by one-to-one matches within two columns. Gradual and uniform families have no designated true boundaries. Report segment counts and lengths. For conserved-region IoU, apply the validation decision threshold to each segment's raw-score mean multiplied by mean nongap coverage over its usable columns. Compare the union of positive segments with the conserved truth on usable coordinates.

## Runtime protocol

Run four wall-time repeats after one warm-up for each measurement. Report median and sample standard deviation. Record Python/NumPy traced allocation peaks in a separate invocation. The peak excludes untracked native workspace and resident input alignment. Record N and L scaling separately, including profile preparation, scoring, contextual processing and the 1,331-point grid. Compare maximum and mean absolute deviations and Spearman on identical inputs. Tests use an absolute tolerance of 10⁻¹².

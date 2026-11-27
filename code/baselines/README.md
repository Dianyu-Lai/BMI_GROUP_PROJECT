# BMI_GROUP_PROJECT

This package implements the earlier entropy, BLOSUM62 and AAindex conservation-scoring experiments for globular protein multiple sequence alignments (MSAs). It requires Python 3.10 or newer and uses the standard library for scoring, benchmarks and tests.

| Path | Contents |
| --- | --- |
| `src/msa_conservation/msa.py` | Alignment input and shared handling of gaps, ambiguous residues, and optional sequence weights. |
| `src/msa_conservation/cli.py` | Command-line entry point for scoring an alignment. |
| `src/msa_conservation/scoring/entropy.py` | Shannon entropy and entropy-based reference scores such as Jensen–Shannon divergence. |
| `src/msa_conservation/scoring/substitution.py` | BLOSUM62-based scoring and substitution-matrix operations. |
| `src/msa_conservation/scoring/physicochemical.py` | Standardized AAindex1 descriptors and probability tilting for BLOSUM62. |
| `src/msa_conservation/scoring/fusion.py` | Combinations of entropy and substitution scores for benchmarking. |
| `src/msa_conservation/scoring/context.py` | Fixed-window and context-adaptive scoring variants for ablation. |
| `../data_generation/synthetic.py` | Rate-driven alignments and frequency-matched within-class and cross-class substitutions. |
| `../benchmarks/metrics.py` | Sensitivity, specificity, and correlation with known conservation patterns. |
| `../benchmarks/run_benchmark.py` | Baseline comparisons, ablations, and runtime and scaling measurements. |
| `data/reference/` | Reference inputs such as BLOSUM62 and selected AAindex properties. |
| `tests/` | Alignment parsing, scoring and synthetic benchmark checks. |

## Run the baseline scores

From `code/baselines` in PowerShell, using Python 3.10 or newer.

```powershell
$env:PYTHONPATH = "..;src"
python -m msa_conservation.cli path\to\alignment.fasta
python -m msa_conservation.cli path\to\alignment.fasta --weights henikoff
python -m benchmarks.run_benchmark --size 64 --length 300 --divergence 0.8 --seed 1
python -m unittest discover -s tests
```

The CLI writes one TSV row per aligned position, numbered from 1. Larger variability scores indicate weaker conservation. Shannon entropy `H` uses a 20-residue alphabet. The BLOSUM62 score `P` is the expected diagonal-normalized substitution cost under two independent draws from the column frequencies. The integrated score is `H * (1 + gamma * M)`, where `M` is the normalized substitution cost conditional on the draws differing. The default `gamma=1` gives the severity term unit weight. The `gamma=0` ablation recovers Shannon entropy.

Frequencies use the 20 canonical residues. A dash or dot is a gap. B, J, O, U, X and Z are treated as noncanonical observations. The TSV reports both nongap coverage and canonical coverage. Sites with zero canonical weight have blank score fields. Uniform weighting is the default. Optional Henikoff weighting reduces the influence of repeated residue types at each position and uses canonical observations to compute sequence weights. The default gap penalty is zero. Setting `--gap-penalty` adds its value times the weighted gap fraction to each of the three variability scores.

The BLOSUM62 reference table contains full-precision, ordered joint residue probabilities. It comes from NCBI BLAST's public-domain `matrix_frequency_data.c`, by Yi-Kuo Yu, Alejandro Schaffer and E. Michael Gertz. The scorer reconstructs half-bit log-odds scores from the row marginals and divides pair costs by the BLOSUM62 probability-weighted mean cost among unlike residue pairs. [NCBI source](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/algo/blast/composition_adjustment/matrix_frequency_data.c#0396) · [Henikoff and Henikoff, 1994](https://pubmed.ncbi.nlm.nih.gov/7966282/)

## Compare reference baselines

```powershell
$env:PYTHONPATH = "..;src"
python -m msa_conservation.cli path\to\alignment.fasta --reference-baselines
python -m benchmarks.run_benchmark --reference-baselines --suite rate
python -m benchmarks.run_benchmark --reference-baselines --suite constraint
```

The flag adds four methods to the existing normalized Shannon entropy, substitution cost, and fusion scores. The optional reference methods use the original BLOSUM62 table at any AAindex strength. Native conservation scores are negated so all columns retain the convention that higher scores indicate more variability.

| Output column | Definition |
| --- | --- |
| `kl_variability` | Negative KL divergence from the column distribution to the BLOSUM62 background, in bits. |
| `jsd_variability` | Negative equal-weight Jensen–Shannon divergence between those distributions, in bits. |
| `blosum62_sop_variability` | Negative weighted mean of the integer BLOSUM62 similarities over distinct sequence pairs. Each pair has weight $w_i w_j$. |
| `adjusted_entropy_variability` | Shannon entropy in nats divided by $2^{(s_{\mathrm{obs}}-s_{\mathrm{bg}})/2}$. The observed mean uses distinct residue types present in the column. The background mean averages their off-diagonal BLOSUM62 row means. |

KL and JSD use canonical frequencies without pseudocounts or window smoothing. Sum-of-pairs omits gaps and ambiguous symbols and gives a blank score when fewer than two sequences supply positive-weight canonical observations. With `--gap-penalty`, the configured additive penalty applies to each available variability score. [Capra and Singh, 2007](https://doi.org/10.1093/bioinformatics/btm270) · [Zhang et al., 2022](https://doi.org/10.1186/s12859-022-04885-7)

## Calibrate BLOSUM62 with AAindex1

```powershell
python -m msa_conservation.cli path\to\alignment.fasta --aaindex-eta 0.5
```

The four descriptors are Kyte–Doolittle hydropathy (`KYTJ820101`), Grantham volume (`GRAR740103`), Klein net charge (`KLEP840101`), and Vihinen mean flexibility (`VINM940101`). Each descriptor is standardized over 20 residues using its population standard deviation. The squared distance `d²(a,b)` is the mean of the four squared descriptor differences. Equal weighting gives each standardized descriptor the same explicit coefficient. The four scales have a maximum absolute Pearson correlation of 0.71. Grantham polarity was left out because its correlation with hydropathy is −0.86 for these 20 residues. [AAindex1 source data](https://www.genome.jp/ftp/db/community/aaindex/aaindex1)

The calibrated joint probability is proportional to `q₀(a,b) exp(-eta d²(a,b)/2)`. Its own row marginals define the half-bit log-odds scores. The cost normalizer remains the original BLOSUM62 value for all strengths, so scores share a common scale. Before clipping, the diagonal-normalized cost increases by `eta d²(a,b) / ln(2)`. The default alignment scorer uses `eta=0`, which gives the original BLOSUM62 scores. The constraint benchmark defaults to `eta=0.5`. A squared distance of 2 multiplies a pair's unnormalized probability by about 0.61 at this strength.

## Run the synthetic benchmarks

The generator draws independent descendants from one ancestral sequence. Three contiguous regions have relative site rates 0.05, 0.25 and 1.0. These settings create separated low, intermediate and high variability regimes. At rate `r` and divergence `d`, each residue changes with probability `1 - exp(-d*r)`. A change draws uniformly from the 19 alternatives. The generator uses relative site rates for Spearman correlation and marks the high-rate region positive for AUROC, sensitivity and specificity. Decision thresholds maximize balanced accuracy on the calibration seed, with ties resolved toward higher thresholds. Test metrics use a separate seed.

The `constraint` suite uses pairs of columns with the same ancestral residue and the same mutation mask. Each pair includes one within-class substitute and one cross-class substitute. The classes are aliphatic `AVLIM`, aromatic `FWY`, uncharged polar `STNQ`, basic `KR` and acidic `DE`. The generator uses these fixed classes directly. Pair members have identical residue counts up to renaming and therefore identical Shannon entropy. The relative mutation rate is 1 throughout, so Spearman in this suite uses the binary cross-class label. Matched-pair accuracy awards half credit for tied pairs. This is a controlled test of class preservation with simplified chemistry and a selected subset of residues.

At a fixed tilt strength and positive fusion coefficient, substitution and fusion scores give the same ordering within each matched pair. Both members have equal entropy and equal unlike-pair probability. AUROC also compares sites across different pairs, where the two scores can order sites differently.

```powershell
python -m benchmarks.run_benchmark --suite constraint --size 64 --length 300 --divergence 0.8 --seed 1
python -m benchmarks.run_benchmark --suite rate --size 64 --length 300 --divergence 0.8 --seed 1 --aaindex-eta 0.5
```

The exploratory sweep used sizes 16 and 64, lengths 90 and 300, divergences 0.3 and 0.8, and calibration seeds 1, 11, 21, 31 and 41. Each test seed was its calibration seed plus 1, giving 40 runs per suite. In the constraint suite, the original BLOSUM62 substitution cost and the AAindex-calibrated version at `eta=0.5` achieved mean AUROCs of 0.842 and 0.905 and matched-pair accuracies of 0.908 and 0.953, respectively. Adjusted entropy achieved a mean AUROC of 0.932 and matched-pair accuracy of 0.921, while Shannon entropy scored 0.500 on both metrics. In the rate suite, Shannon entropy had the highest mean Spearman correlation (0.892). These values describe the controlled synthetic generators.

For `N` sequences, `L` sites and alphabet size `A=20`, computing column frequencies takes `O(NL)` and substitution scoring takes at most `O(LA²)`. Substitution scoring visits only pairs of observed residue types. The CLI streams scores with `O(NL + A²)` memory including the alignment and reference matrix. The benchmark stores its scores and takes `O(L log L)` additional time for ranking and threshold selection.

AAindex calibration builds one global matrix per strength, with `O(A²d)` setup for `d` descriptors. The reference sum-of-pairs groups pairs by residue type and uses a second pass over each column to accumulate within-type weights. Its total time is `O(NL + LA²)`. KL and JSD add `O(LA)`, and adjusted entropy adds up to `O(LA²)`.

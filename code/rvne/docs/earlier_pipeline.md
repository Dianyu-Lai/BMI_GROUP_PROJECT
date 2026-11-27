# BMI_GROUP_PROJECT

BMI3 Algorithms ICA Group 1 project on position-specific conservation scoring for globular protein multiple sequence alignments (MSAs). The core scoring methods will be implemented in Python 3.

| Path | Planned responsibility |
| --- | --- |
| `src/msa_conservation/msa.py` | Alignment input and shared handling of gaps, ambiguous residues, and optional sequence weights. |
| `src/msa_conservation/cli.py` | Command-line entry point for scoring an alignment. |
| `src/msa_conservation/scoring/entropy.py` | Shannon entropy and entropy-based reference scores such as Jensen–Shannon divergence. |
| `src/msa_conservation/scoring/substitution.py` | BLOSUM62-based scoring and substitution-matrix operations. |
| `src/msa_conservation/scoring/physicochemical.py` | Standardized AAindex1 descriptors and probability tilting for BLOSUM62. |
| `src/msa_conservation/scoring/fusion.py` | Combinations of entropy and substitution scores for benchmarking. |
| `src/msa_conservation/scoring/context.py` | Fixed-window and context-adaptive scoring variants for ablation. |
| `benchmarks/synthetic.py` | Rate-driven alignments and frequency-matched within-class and cross-class substitutions. |
| `benchmarks/metrics.py` | Sensitivity, specificity, and correlation with known conservation patterns. |
| `benchmarks/run_benchmark.py` | Baseline comparisons, ablations, and runtime and scaling measurements. |
| `data/reference/` | Reference inputs such as BLOSUM62 and selected AAindex properties. |
| `data/real/` | Optional validation alignments from globular protein families. |
| `data/synthetic/` | Generated alignments and their ground-truth labels. |
| `results/` | Generated benchmark tables and figures. |
| `tests/` | Future tests and fixtures for alignment parsing, scoring, and benchmark inputs. |
| `docs/` | Method, benchmark design, and complexity notes. |

## Run the baseline scores

From the repository root in PowerShell, using Python 3.10 or newer:

```powershell
$env:PYTHONPATH = "src"
python -m msa_conservation.cli path\to\alignment.fasta
python -m msa_conservation.cli path\to\alignment.fasta --weights henikoff
python -m benchmarks.run_benchmark --size 64 --length 300 --divergence 0.8 --seed 1
python -m unittest discover -s tests
```

The CLI writes one TSV row per aligned position, numbered from 1. Larger variability scores indicate weaker conservation. Shannon entropy `H` uses a 20-residue alphabet. The BLOSUM62 score `P` is the expected diagonal-normalized substitution cost under two independent draws from the column frequencies. The integrated score is `H * (1 + gamma * M)`, where `M` is the normalized substitution cost conditional on the draws differing. The default `gamma=1` gives the severity term unit weight. The `gamma=0` ablation recovers Shannon entropy.

Frequencies use the 20 canonical residues. A dash or dot is a gap. B, J, O, U, X and Z are treated as noncanonical observations. The TSV reports both nongap coverage and canonical coverage. Sites with zero canonical weight have blank score fields. Uniform weighting is the default. Optional Henikoff weighting reduces the influence of repeated residue types at each position and uses canonical observations to compute sequence weights. The default gap penalty is zero. Setting `--gap-penalty` adds its value times the weighted gap fraction to each of the three variability scores.

The BLOSUM62 reference table contains full-precision, ordered joint residue probabilities. It comes from NCBI BLAST's public-domain `matrix_frequency_data.c`, by Yi-Kuo Yu, Alejandro Schaffer and E. Michael Gertz. The scorer reconstructs half-bit log-odds scores from the row marginals and divides pair costs by the BLOSUM62 probability-weighted mean cost among unlike residue pairs. [NCBI source](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/algo/blast/composition_adjustment/matrix_frequency_data.c#0396) · [Henikoff and Henikoff, 1994](https://pubmed.ncbi.nlm.nih.gov/7966282/)

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

The exploratory sweep used sizes 16 and 64, lengths 90 and 300, divergences 0.3 and 0.8, and calibration seeds 1, 11, 21, 31 and 41. Each test seed was its calibration seed plus 1. At `eta=0.5`, the mean constraint AUROC was 0.842 for BLOSUM62 and 0.905 for calibrated BLOSUM62. Mean matched-pair accuracy rose from 0.908 to 0.953. Shannon entropy retained the highest mean Spearman correlation in the rate suite (0.892). The strength was held fixed within each sweep, and these results describe the controlled generator.

For `N` sequences, `L` sites and alphabet size `A=20`, computing column frequencies takes `O(NL)` and substitution scoring takes at most `O(LA²)`. Substitution scoring visits only pairs of observed residue types. The CLI streams scores with `O(NL + A²)` memory including the alignment and reference matrix. The benchmark stores its scores and takes `O(L log L)` additional time for ranking and threshold selection.

AAindex calibration builds one global matrix per strength, with `O(A²d)` setup for `d` descriptors. Context adaptation and published-method comparisons are planned extensions.

# SCAW-RVNE v2: penalized segmentation with fixed local context

Version `scaw-v2` fits changepoints and applies fixed-radius windows within each segment. See the [current benchmark](../results/scaw_revision/README.md) for evaluation.

## Objective and one context parameter

Let y_i be raw base-20 RVNE and w_i the nongap fraction at a scorable column; otherwise w_i=0. Coverage is a reliability proxy, not a sequence count or an effective sample size. Define

\[
\mu_{a:b}=\frac{\sum_{i=a}^{b-1}w_i y_i}{\sum_{i=a}^{b-1}w_i},\qquad
C(a,b)=\sum_{i=a}^{b-1}w_i(y_i-\mu_{a:b})^2.
\]

The implemented objective is

\[
\min_{\tau}\left\{\sum_{j=0}^{K}C(\tau_j,\tau_{j+1})+\beta K\right\},
\quad\beta=c\,\widehat\sigma^2\,2\log(\max(2,m)).
\]

Here K counts changepoints, m counts scorable columns, and **c=`penalty_scale` is the only adjustable context parameter**, default 1. The Gaussian working model is Var(y_i)=sigma²/w_i with constant regime mean. Under this model C/sigma² is the variable part of minus twice the log likelihood. Counting an extra mean and a boundary gives the chosen 2 log(m) BIC-style scale. Treating a discrete searched boundary as an ordinary parameter is an approximation: this is **our BIC-style plug-in adaptation**, not the original BIC theorem, modified BIC, a calibrated significance level, or a formula validated specifically for RVNE.

On the set A of pairs of adjacent original columns that are both scorable,

\[
\widehat\sigma^2=\frac1{|A|}\sum_{i\in A}
\frac{(y_{i+1}-y_i)^2}{1/w_i+1/w_{i+1}}.
\]

With all w=1 and no missing sites this is Rice's first-difference estimator. The weighted form is our derivation under the working variance model. Changes add positive bias; positive residual serial correlation can add negative bias. Missing runs are never compressed when estimating noise. If A is empty, use full-profile coverage-weighted variance; a single usable point has variance zero. This fallback is deterministic and is not claimed as a noise estimator with a consistency guarantee. No fitted noise floor is added.

## Exact solver and missing data

Centered weighted prefix sums give C(a,b) in O(1). PELT uses F(0)=-beta and

\[
F(t)=\min_{s<t}\{F(s)+C(s,t)+\beta\}.
\]

Weighted SSE satisfies C(a,t) >= C(a,s)+C(s,t), so the pruning constant is zero. Optimize on scorable columns only, allow at least one scorable column per segment, and impose no maximum. Missing columns provide no split evidence: map a changepoint to the first scorable column of the right segment; interior missing runs belong to the left segment, leading missing columns to the first. Boundary positions inside missing runs are unidentifiable. Segments cover the original coordinates and have positive support. A completely unscorable profile has no segments. Missing center scores remain NaN in every method.

Prefer fewer segments on numerical ties, then the earliest predecessor. The tolerance is 64 machine epsilons times the profile centered sum of squares or penalty, whichever is larger; it only handles prefix-sum roundoff. Exhaustive partition tests check optimal objective values on randomized short weighted profiles, including zero penalty. PELT has worst-case O(m²) time and O(L+m) space. Its expected linear cost requires additional distribution/change-frequency conditions and is not guaranteed for MSA profiles.

## Fixed window

Set radius=3 and center weight=1/2 as fixed design constants informed by Capra & Singh (2007), equation (6). For a scorable center i, let N_i be the scorable neighbors, excluding i, within [i-3,i+3] and its assigned segment:

\[
\widehat y_i=\tfrac12 y_i+\tfrac12
\frac{\sum_{j\in N_i}w_jy_j}{\sum_{j\in N_i}w_j},\quad
S_i=w_i\widehat y_i.
\]

When there are no usable neighbors, use y_i. Coverage weighting, explicit neighbor exclusion, missing-data handling, and segment clipping are our implementation rules; this is not an exact reproduction of the paper's full conservation pipeline. The paper does not establish radius 3 as optimal for RVNE. Clipping preserves short motifs **only if the segmentation detects them**. There is no lambda objective, candidate-radius search, support cutoff, or adjustable min/max segment length.

`window_size` is the coordinate span of the clipped window, including missing coordinates. `window_support` is its sum of usable coverage, including the center. `local_variance` is the coverage-weighted variance in that window around its weighted mean, not the sampling variance of the blended estimator. The segmentation ablation emits full-segment means and support. Uniform and Gaussian ablations use fixed radius 3; the Gaussian sigma is 1.5. Their window constants are design choices, not new UI parameters.

## Interfaces and reproducibility

CLI: `--penalty-scale 1`; JSON: `"context": {"penalty_scale": 1.0}`. Settings downloads contain `context_version`, the effective noise variance and penalty, fixed radius and center weight, and solver identity. `context_version` records the formula contract. Legacy settings and flags fail explicitly; no numerical conversion from old beta is valid because the scaling changed.

The frontend uses `results/scaw_revision/parameters.json` for its study preset and reads the benchmark summary from the same directory. The evaluation uses fixed RVNE values from the original training and validation split.

Limitations: the MSA score profile is neither demonstrably independent Gaussian data nor guaranteed inverse-coverage-variance data. Phylogeny, residue-specific background effects, alignment errors, gradients and heavy gaps can defeat these assumptions. These are estimated score regimes, not motifs with biological significance or causal annotations. Fixed binary thresholds learned on historical simulations did not transfer uniformly to the fresh INDELible set. The benchmark report gives the results by cohort.

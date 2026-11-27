# Segmentation-Constrained Adaptive-Window RVNE

SCAW-RVNE estimates conservation at every usable alignment column from the original raw RVNE profile D. Let w_i be the nongap fraction for a usable site and zero for all missing or excluded sites. A prefix of w, w(D−m), and w(D−m)² gives the support C, weighted mean μ, and within-interval weighted SSE in constant time. Here m is the profile's weighted mean. Centering limits cancellation. For roundoff, the prefix-based SSE is bounded below by zero.

## Segmentation

Minimize the sum of the coverage-weighted SSE plus βV per segment, where V is the variance on the full usable profile. Scale normalization makes β comparable across RVNE settings. The recurrence for the best prefix cost is F(t) = min over feasible s of [F(s) + SSE(s,t) + βV], with F(0) = 0. Backtrace the minimizing split points. At exact ties prefer fewer segments, then an earlier start. Segments require positive valid support. Set the minimum length to three by default to limit single-site fragmentation, using the whole alignment when it is shorter. There is no default maximum. A user-supplied maximum can save computation at the cost of introducing splits in long homogeneous regions. Gappy positions have zero segmentation weight, so a boundary inside a gap run may be weakly identifiable. Missing scores remain missing.

## Windows and integration

For each scorable column i enumerate centered radii 0, 1, 2, 4 and 8. Clip each candidate to the alignment and the assigned segment. Select the candidate that minimizes local variance / V + λ / C. Prefer larger support on exact objective ties. A constant profile has zero variance ratio. The minimum candidate support is 1.5 by default. Radius zero remains available when the center has low coverage or belongs to a very short segment. The variance term reduces smoothing over heterogeneous neighborhoods. The reciprocal support term rewards wider windows in homogeneous regimes. Candidate limits and λ are ablated because radius-zero bias and oversmoothing are both possible. A local estimate remains missing if its center was excluded. The final score is the local weighted mean multiplied by the center's coverage, keeping the same center penalty as original RVNE.

    D, c, valid <- score original RVNE
    prefixes <- weighted moments of valid (D, c)
    segments <- optimal_partition(prefixes, beta, min_length, max_length)
    for each usable column i:
        windows <- candidates clipped to i's segment
        estimate, width, variance, support <- minimize interval objective
        final_i <- c_i * estimate
    export site and segment tables

The ablations reuse the same raw scores and gap mask. Fixed uniform and Gaussian windows use the same maximum radius (default 8). Gaussian uses σ = radius/2. The segmentation-only ablation uses the segment mean. The adaptive-only ablation permits windows throughout the alignment. Segment outputs use the raw profile, with optional user-provided lower and upper raw-score thresholds for variable, intermediate and conserved classes. The default is unclassified.

## Reference computation and tuning

Keep the equation-11 reference path. The equivalent trace path caches log σ, computes observed eigenvalues in bounded batches, and contracts ρ against log σ. For a parameter grid, cache observed eigenvalue entropies per s1 within each column batch and cached log σ per (r,s2). Score with no calibration data through the command line. In the experiments, select RVNE using the published grid plus a local refinement using training families, select β and λ using separate validation families, and lock test predictions and classification cutoffs before evaluating test truth. The benchmark protocol records the exact selection rules.

## Prior methods and limitations

Capra and Singh (2007) incorporated a fixed local context window into protein conservation prediction. Caffrey et al. (2004) used substitution similarity in VNE. Adaptive weights smoothing selects local homogeneous neighborhoods, bilateral filtering downweights cross-edge neighbors, and Killick et al. describe exact penalized changepoint detection. These are established ideas. The contribution investigated here is the use of coverage-weighted changepoints as hard constraints on independently selected RVNE window sizes, with ablations for both components. The reviewed papers did not describe this exact combination. The search is limited and does not establish a priority claim.

Gaussian and uniform averaging can blur sharp edges. SCAW-RVNE can instead lose short motifs when segmentation merges them or oversegment random fluctuations at low depth. Smooth gradients can be represented as steps. Heavy gaps weaken support and can hide true changes. RVNE is sensitive to the identity of a conserved residue through the background, so high accuracy on mutation-rate truth is an empirical question. The segment classes represent estimated conservation regimes.

- [Capra and Singh, 2007](https://doi.org/10.1093/bioinformatics/btm270)
- [Caffrey et al., 2004](https://doi.org/10.1110/ps.03323604)
- [Polzehl and Spokoiny, 2000](https://doi.org/10.1111/1467-9868.00235)
- [Tomasi and Manduchi, 1998](https://doi.org/10.1109/ICCV.1998.710815)
- [Killick, Fearnhead and Eckley, 2012](https://doi.org/10.1080/01621459.2012.737745)

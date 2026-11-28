# Time and memory complexity

Let N be sequence count, L alignment length, K alphabet size, M usable columns, R fixed window radius (3), B batch size, and G = G_r G_1 G_2 parameter settings. Memory excludes the underlying NumPy/LAPACK runtime and library caches unless stated.

| Stage | Time | Additional space |
| --- | --- | --- |
| FASTA validation and input | O(NL) | O(NL) |
| Henikoff weights, profiles and coverage | O(NL + LK) | O(N + LK) |
| One A(s) with PSD check | O(K³) | O(K²) |
| One density construction | O(K²) | O(K²) |
| Reference equation 11 over M sites | O(MK³) | O(L + K²) beyond profiles |
| Cached log σ setup | O(K³) | O(K²) per cached background |
| Cached trace path | O(MK³ + K³) | O(L + BK²) beyond profiles |
| Grid with reuse | O((G₁M + GᵣG₂)K³ + GMK²) | O(GM + (B + GᵣG₂)K²) beyond profiles |
| Coverage and score prefix sums | O(L) | O(L) |
| PELT weighted SSE partition and backtrace | O(M²) worst case; expected O(M) only under PELT conditions | O(L + M) |
| Fixed boundary-clipped context | O(L) | O(L) |
| Fixed uniform and Gaussian windows | O(L), O(LR) | O(L), O(L + R) |
| Site ranking metrics | O(M log M) | O(M) |
| Matching sorted boundary sets | O(P + T) for P predicted and T true boundaries | O(1) |

The reference builds one background density and recomputes its eigendecomposition inside each independent equation-11 call. Observed eigenvectors and dense overlap also cost O(K³). The optimized path removes the repeated background decompositions and observed eigenvectors. Its symmetric eigenvalue computation remains cubic in K. Batching trades traced peak memory for a constant-factor time reduction.

Running equation 11 independently for each grid setting costs O(GMK³). Reuse reduces eigendecompositions to O(G₁M + GᵣG₂), while every setting still pays for its dense cross term. All G×M scores are returned for rank metrics, so peak grid memory grows with L even though eigensolver work is performed in batches.

SCAW v2 computes noise from adjacent usable original columns in O(L), then solves weighted SSE + beta K by PELT on the M usable columns. Prefix statistics cost O(L), and each interval query costs O(1). PELT removes predecessors using weighted SSE superadditivity. Long homogeneous regions may still require quadratic work; the expected-linear theorem is conditional and is not a general guarantee for MSA profiles. There is no maximum-length cap or radius search.

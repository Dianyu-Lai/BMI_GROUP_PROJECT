# Original RVNE specification

Johansson and Toh (2010), [Relative von Neumann entropy for evaluating amino acid conservation](https://doi.org/10.1142/S021972001000494X), equations 8–11 and sections 2–3.

## Input and weights

The input is a rectangular protein alignment in FASTA format. The canonical alphabet is `ARNDCQEGHILKMFPSTWYV`. Henikoff position-based weights sum 1/(d n_a) over each canonical residue observation, where d is the number of residue types and n_a is that residue's count in the column. Normalize the sequence weights globally, then normalize weighted canonical counts within each column. An alignment without any canonical observations has no scorable sites. The parser treats dash and dot as alignment gaps. B, J, O, U, X and Z count as nongap observations for coverage and contribute no mass to p or to sequence weights. The paper does not state an ambiguity policy, so this last rule is an implementation convention.

## Matrices and scores

Read the original four-decimal BLOSUM62 joint pair frequencies from [the AUGUSTUS distribution](https://github.com/Gaius-Augustus/Augustus/blob/master/config/profile/blosum62.qij), Git blob e3e09d03009c54e22903f85de250cff1d8284d12. The file stores a lower triangle in the above amino acid order. Mirror it to form M. Its rounded entries sum to 0.9987, so divide the row totals by 0.9987 for the normalized background. This rescaling cancels in C_ab = M_ab / sqrt(M_aa M_bb). Conditioning the table on L, M and W reproduces the displayed precision of Johansson and Toh's Table 1. The separate existing high-precision NCBI table remains in use for the substitution-severity benchmark.

Use A(0) = I, A(1) = 11ᵀ and A(s)_ab = C_ab^((1−s)/s) inside those endpoints. Use q = 1/20 + r(q̄ − 1/20), ρ = diag(sqrt(p)) A(s1) diag(sqrt(p)), and σ = diag(sqrt(q)) A(s2) diag(sqrt(q)). Diagonals of A are one, so both density matrices have trace one. The original grid has 11 values per parameter. r and s1 span 0 to 1 in tenths. s2 spans 0 to 0.9 in tenths, plus 0.99. At s2 = 1 the background has rank one and its logarithm is singular. Arbitrary observed densities need not share its support.

When (λ_i,v_i) and (μ_j,w_j) are eigenpairs of ρ and σ, compute equation 11 as

    D = Σ_{λ_i>0} λ_i [log_20 λ_i − Σ_j |v_iᵀw_j|² log_20 μ_j].

The independent trace calculation is Σ_i λ_i log_20 λ_i − tr(ρ log_20 σ). Use zero for terms with λ_i = 0. Reject materially negative eigenvalues and nonpositive background eigenvalues. Roundoff-negative observed eigenvalues within 1e−12 are treated as numerical zeros in the scalar entropy term. No matrix is projected onto the PSD cone. An exact uniform background with s2=0 reduces the score to 1−VNE. With s1 also zero it becomes 1−Shannon entropy. With r>0 and s1=s2=0 it becomes relative Shannon entropy.

Preserve all alignment coordinates. For gap fraction g/N > 0.30 or zero canonical observation mass, the score is missing. Otherwise report raw divergence D and paper score (1−g/N)D. Coverage is the unweighted fraction of nongap sequence rows. All outputs orient higher values toward greater conservation. Relative entropy can exceed one for nonuniform or correlated backgrounds. Defaults are (r,s1,s2) = (0.8,0.8,0.8).

## Reproduction scope

The implementation verifies the equations, Table 1, parameter grid, cutoff and special cases. Reproducing the published catalytic-site AUC would require the original alignment and annotation snapshot.

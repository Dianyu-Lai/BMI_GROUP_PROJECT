# MAFFT and ClustalW for real-data validation

Use MAFFT L-INS-i as the primary aligner for a small set of related globular protein domains. Use ClustalW 2.1 as a progressive alignment control. Run both on the same sequence set with the same RVNE or SCAW-RVNE settings. This experiment tests how sensitive conservation scores are to alignment choice. The synthetic benchmark supplies known conservation and variable-region labels for evaluating the scoring method.

| Aligner or mode | Method and role |
| --- | --- |
| MAFFT `linsi` | Combines local pairwise alignments with iterative refinement. Use for a small family with a shared domain and possible insertions. The official site suggests up to roughly 200 sequences and 2,000 sites for this mode. |
| MAFFT `ginsi` | Uses global pairwise alignments with iterative refinement. Use for similar-length sequences that are homologous end to end. |
| MAFFT `einsi` | Uses local pairwise alignment with a generalized affine gap cost and iterative refinement. Use when homologous blocks are separated by long insertions or unalignable regions. |
| MAFFT `auto` or `fftns2` | `auto` selects a strategy by input size. `fftns2` fixes the fast progressive strategy for a reproducible runtime comparison with larger inputs. |
| ClustalW 2.1 | Builds a guide tree and aligns sequences progressively, using sequence weights, position-specific gap penalties, and substitution matrix choices. Early guide-tree merges can affect later columns. Use the default alignment settings as a historical control. |

The expensive pairwise stage in the accuracy-oriented MAFFT modes performs up to $\binom{N}{2}$ dynamic-programming alignments. For sequences of similar length $L$, that stage is roughly $O(N^2L^2)$ in time before progressive alignment and refinement. ClustalW also performs pairwise comparisons to build its guide tree, then merges sequence groups. Report each aligner's observed time and dimensions with the strategy and thread count. The project runs both aligners externally and reuses the project's Python conservation scorer. The ClustalW consensus symbols in the synthetic benchmark are a separate column-scoring baseline.

## Run the example

Install [MAFFT 7](https://mafft.ddbj.nig.ac.jp/alignment/software/) and [ClustalW 2.1](https://www.theclustal.com/clustal2/). On Windows, use the [MAFFT Windows package](https://mafft.ddbj.nig.ac.jp/alignment/software/windows_without_cygwin.html) and the ClustalW command-line executable. If the executables are outside `PATH`, pass their full paths.

Run from `code/rvne` after installing `requirements.txt`.

```powershell
$env:PYTHONPATH = "src"
python -m msa_conservation.cli real-data --input data/real/globins.fasta --reference-id P02144 --mafft-mode linsi --output-dir results/globins
```

Use `--mafft-exe 'C:\path\mafft-win\mafft.bat'` and `--clustalw-exe 'C:\path\clustalw2.exe'` for explicit paths. Use `--aligner mafft` or `--aligner clustalw` for a single aligner. Use `--method rvne` for the original RVNE or `--settings results/scaw_revision/parameters.json` to reuse the frozen RVNE and SCAW parameters. For a matched alignment sensitivity check, compare the same scoring method and settings across aligners. Use a new output directory for each run.

The input is an unaligned FASTA with at least two nonempty sequences. Headers must have unique first tokens. Use the 20 canonical amino acids. Gap characters, ambiguous residues and stop symbols are rejected explicitly so both tools receive the same sequences. For a larger family, select homologous globular domains with comparable boundaries, document exclusions and redundancy filtering, and retain the accession and start coordinate of each domain. The reference defaults to the first sequence. Coordinates start at 1 in the sequence supplied, so a domain sliced from a full-length entry requires the source start offset for annotation lookup.

## Read the outputs

| File | Contents |
| --- | --- |
| `<aligner>.fasta` | Alignment in input sequence order, with original FASTA IDs restored. |
| `<aligner>.scores.tsv` | Existing RVNE score columns for every alignment column, including columns with gaps in the reference. |
| `<aligner>.segments.tsv` | SCAW segments in that aligner's column coordinates. Available with `scaw-rvne`. |
| `reference_scores.tsv` | One row per residue of the ungapped reference sequence, with each aligner's original column, nongap coverage and score. Blank scores retain the score filter. |
| `summary.json` | Input SHA-256, sequence count, reported aligner versions, commands, alignment wall times, lengths, unweighted mean gap fraction and numbers of scored reference residues. Dual-aligner runs also report Pearson correlation across jointly scored reference positions. Pearson is null when fewer than two residues are scored by both or either compared profile is constant. |
| `scoring_settings.json`, `<aligner>.log` | Effective scorer settings and aligner diagnostics. |

`data/real/globins.fasta` contains five canonical UniProt sequences, with only headers and line wrapping changed. The entries are human myoglobin [P02144](https://www.uniprot.org/uniprotkb/P02144/entry), horse myoglobin [P68082](https://www.uniprot.org/uniprotkb/P68082/entry), sperm whale myoglobin [P02185](https://www.uniprot.org/uniprotkb/P02185/entry), human hemoglobin alpha [P69905](https://www.uniprot.org/uniprotkb/P69905/entry), and human hemoglobin beta [P68871](https://www.uniprot.org/uniprotkb/P68871/entry). The aligners process the same five sequences without gaps supplied. In the default reference sequence P02144, UniProt annotates His94 as the proximal haem iron ligand. Inspect residue 94 in `reference_scores.tsv` along with neighbouring sites and coverage.

The five globins form a worked example. For a quantitative alignment accuracy test, use a versioned structure-informed reference such as BAliBASE 3.0 and compare aligned residue pairs and exact columns in its annotated core blocks. For a broader scoring validation, select several curated globular-domain alignments from Pfam and map independently annotated residues to the reference sequence. Report coverage and the ranks of annotated sites. Unannotated residues may also be functionally important, and the BAliBASE core blocks encode alignment reliability. Use the synthetic test labels to quantify sensitivity and specificity for known variable regions.

## Sources

- Katoh and Standley (2013), [MAFFT v7](https://doi.org/10.1093/molbev/mst010) and the [MAFFT manual](https://mafft.ddbj.nig.ac.jp/alignment/software/manual/manual.html).
- Thompson, Higgins and Gibson (1994), [ClustalW](https://doi.org/10.1093/nar/22.22.4673) and Larkin et al. (2007), [Clustal W and X 2.0](https://doi.org/10.1093/bioinformatics/btm404).
- Thompson et al. (2005), [BAliBASE 3.0](https://pubmed.ncbi.nlm.nih.gov/16044462/) and [Pfam globin family PF00042](https://www.ebi.ac.uk/interpro/entry/pfam/PF00042/).

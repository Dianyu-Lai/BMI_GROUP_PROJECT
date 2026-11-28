# Published conservation baselines

This package scores columns of an aligned protein FASTA using Python 3.10 or newer and the standard library.

| Path | Purpose |
| --- | --- |
| `src/msa_conservation/msa.py` | FASTA parsing, gap and ambiguous-residue handling, sequence weighting and column frequencies. |
| `src/msa_conservation/cli.py` | Command-line entry point for scoring an alignment. |
| `src/msa_conservation/scoring/baseline.py` | Per-column scoring and coverage reporting. |
| `src/msa_conservation/scoring/entropy.py` | Shannon entropy, KL divergence, Jensen–Shannon divergence and adjusted entropy. |
| `src/msa_conservation/scoring/substitution.py` | BLOSUM62 background frequencies, integer scores and weighted sum-of-pairs scoring. |
| `data/reference/blosum62_joint.tsv` | BLOSUM62 joint probabilities from [NCBI BLAST](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/algo/blast/composition_adjustment/matrix_frequency_data.c#0396). |
| `tests/` | Parsing, numerical checks and benchmark integration tests. |

## Score an alignment

Run from `code/baselines` in PowerShell. In a POSIX shell, use `export PYTHONPATH=..:src`.

```powershell
$env:PYTHONPATH = "..;src"
python -m msa_conservation.cli ../rvne/data/examples/small_alignment.fasta
python -m msa_conservation.cli ../rvne/data/examples/small_alignment.fasta --weights henikoff --reference-baselines
python -m unittest discover -s tests
```

The CLI writes one TSV row per aligned position, numbered from 1. Larger scores indicate more variability. The default reports normalized Shannon entropy. The `--reference-baselines` flag adds four scores.

| TSV score | Definition |
| --- | --- |
| `entropy_variability` | Shannon entropy divided by the natural log of 20. |
| `kl_variability` | Negative KL divergence from the column frequencies to the BLOSUM62 background, in bits. |
| `jsd_variability` | Negative equal-weight Jensen–Shannon divergence between the column frequencies and BLOSUM62 background, in bits. |
| `blosum62_sop_variability` | Negative mean integer BLOSUM62 similarity over distinct sequence pairs, weighted by the product of their sequence weights. |
| `adjusted_entropy_variability` | Shannon entropy in nats divided by a BLOSUM62 similarity coefficient computed from the distinct residue types present in the column. |

Frequencies are normalized over the 20 canonical amino acids. A dash or dot is a gap. B, J, O, U, X and Z count toward nongap coverage but are excluded from the frequency distribution. The CLI reports both nongap and canonical coverage. All score fields are blank when canonical weight is zero. The sum-of-pairs field also requires at least two positive-weight canonical observations.

Sequences are equally weighted by default. Select `--weights henikoff` to downweight redundant residues. The optional `--gap-penalty` adds its value times the weighted gap fraction to each available score and defaults to zero.

## Evaluate on synthetic alignments

From the same directory, run the shared baseline benchmark with the desired synthetic suite.

The rate suite uses NumPy and [INDELible 1.03](../data_generation/README.md). Set `--data-output` to save its alignments, truth and raw run files. The matched-constraint suite uses the standard library to construct its within-class and cross-class controls.

```powershell
python -m benchmarks.run_benchmark --suite rate --reference-baselines
python -m benchmarks.run_benchmark --suite constraint --reference-baselines
```

The rate suite compares rankings across low-, medium- and high-rate regions. The constraint suite matches mutation masks and frequencies between within-class and cross-class substitution columns. Both fit decision thresholds on a calibration alignment and report metrics on a separate test alignment. See the [benchmark guide](../benchmarks/README.md) for RVNE and SCAW-RVNE comparisons.

For `N` sequences, `L` columns and an alphabet of size `A=20`, column profiling takes `O(L(N+A))` time. Shannon, KL and JSD take `O(A)` per column. BLOSUM62 sum-of-pairs groups pairs by residue type and uses at most `O(N+A^2)` time per column, including the within-type weight pass. Adjusted entropy takes at most `O(A^2)` per column.

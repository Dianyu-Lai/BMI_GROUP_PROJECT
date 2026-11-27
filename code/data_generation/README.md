# Synthetic protein alignments

Run from `code/` in PowerShell with Python 3.10 or newer. Install NumPy from [the RVNE requirements](../rvne/requirements.txt). The generators use the scorers' shared alignment format.

```powershell
python -m pip install -r rvne/requirements.txt
$env:PYTHONPATH = "rvne/src"
```

In a POSIX shell, use `export PYTHONPATH=rvne/src`.

| Path | Contents |
| --- | --- |
| `rvne_simulation.py` | Reference RVNE scenarios and single-family FASTA and site-truth export. |
| `generate_pattern_data.py` | Supplemental A–D dataset suites and manifests. |
| `synthetic.py` | Rate-driven and matched-pair synthetic alignments for the earlier baseline comparisons. This generator uses the standard library. |
| `datasets/` | Local generated or extracted inputs. |

## Generate a single family

```powershell
python -m data_generation.rvne_simulation --scenario motifs --seed 42 --output data_generation/datasets/motifs.fasta --truth data_generation/datasets/motifs_truth.tsv
```

## Generate supplemental test datasets

If you already have the supplemental ZIP archives, extract them under `datasets/` so that `patterns_abc/manifest.tsv` and `patterns_d/manifest.tsv` sit alongside their FASTA and truth files. The commands below regenerate the inputs from their recorded seeds.

```powershell
python -m data_generation.generate_pattern_data --output data_generation/datasets/patterns_abc
```

Choose a new output directory. The generator writes 36 aligned FASTA files, matching site-truth TSVs, and a manifest with the settings and seeds. Each alignment has 50 sequences and 240 columns. Each setting has four independent replicates.

| Pattern | Site-rate profile | Varied setting |
| --- | --- | --- |
| A | Uniformly conserved at rate 0.06 | Divergence 0.35, 0.7, 1.4 |
| B | Conserved flanks at 0.06 and a central variable block at 1.0 | Variable block fraction 25%, 40%, 60% |
| C | A central conserved motif at 0.06 in a variable background at 1.0 | Motif length 3, 5, 9 residues |

The truth TSV records the one-based position, simulated mutation rate, expected conservation probability, conserved label and a flag for a boundary immediately before the position. The supplemental seeds 50100–50135 are independent of the reference benchmark's seeds.

Generate the optional three-rate Pattern D in a new directory.

```powershell
python -m data_generation.generate_pattern_data --suite d --output data_generation/datasets/patterns_d
```

Pattern D uses five consecutive blocks with low, medium, high, medium and low mutation rates (0.06, 0.32 and 1.0). The central high-rate block takes up 25%, 40% or 60% of the alignment, matching Pattern B's settings. The medium-rate blocks take up half of each flank, rounded down to a whole number of columns. At divergence 0.7, the three per-descendant mutation probabilities are approximately 4.1%, 20.1% and 50.3%, using `1 - exp(-divergence * rate)`. The suite contains 12 families (four seeds per setting, 50136–50147), with 50 sequences and 240 columns each.

Use Spearman correlation against the `conservation` column to check whether scores recover the three-tier ordering. The binary `conserved` column marks only the low-rate sites as conserved. The truth file records all four rate transitions for boundary evaluation.

### Pattern D with rate-specific gaps

```powershell
python -m data_generation.generate_pattern_data --suite d-gaps --output data_generation/datasets/patterns_d_gaps
```

| Rate regime | Substitution rate | Gap probability per aligned cell |
| --- | ---: | ---: |
| Low | 0.06 | 3% |
| Medium | 0.32 | 12% |
| High | 1.00 | 25% |

Each of the 12 families matches a gap-free D family in substitution seed, aligned positions and substitution-rate truth. After generating substitutions, the simulator independently masks each sequence and column as `-` according to its rate regime. This models position-dependent missingness in an existing alignment. The 25% setting places the high-rate block near the default 30% scoring threshold, allowing some finite-sample columns to exceed it. The manifest records both seeds and all three gap probabilities. Each truth TSV adds the expected gap probability and observed gap fraction for each column.

Use the `patterns_d_gaps` manifest in a separate benchmark run when comparing against gap-free D. Both manifests use the same pattern label so that the benchmark applies the D-specific diagnostics to both.

See the [benchmark commands](../benchmarks/README.md) to evaluate generated inputs.

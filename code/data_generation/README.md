# Synthetic protein alignments

INDELible 1.03 evolves globular-protein sequences along a rooted tree under an empirical amino-acid substitution model. The simulator tracks the history of insertions and deletions and exports the true alignment.

## Install and configure

Install Python 3.10 or newer, NumPy from [`rvne/requirements.txt`](../rvne/requirements.txt), and **INDELible 1.03** as an external program. Windows users can extract the Windows executable from the [INDELible distribution](http://abacus.gene.ucl.ac.uk/software/indelible/) or build its C++ source with MinGW and `-DWINDOWS`. On Debian or Ubuntu, `sudo apt install indelible` installs the packaged program where available. Other Linux platforms can compile the release source using its build instructions.

The generator checks `--indelible-bin` first and `INDELIBLE_BIN` second. With neither set, it looks in `.tools/indelible/` under the repository root for `indelible.exe` on Windows or `indelible` on Linux, then checks `PATH`. The project-local binary directory is Git-ignored, so a fresh checkout needs its own executable. From the repository root, `python code/data_generation --check` verifies the Python imports and runs a tiny simulation. See the [assessor walkthrough](../frontend/README.md#generate-and-score-your-own-test-alignment) for the full generation and scoring example.

The `python -m` commands below use the project's `code/` directory. Run the following commands from that directory in PowerShell. For a POSIX shell use `export PYTHONPATH=rvne/src` and `export INDELIBLE_BIN=/path/to/indelible` when needed.

```powershell
python -m pip install -r rvne/requirements.txt
$env:PYTHONPATH = "rvne/src"
$env:INDELIBLE_BIN = "C:\tools\INDELible\indelible.exe"
```

## Model and truth

Each family shares one seeded, rooted topology and taxon set across its regions. The base tree's root-to-tip distance is `--divergence` in **expected substitutions per site** for a region with rate 1. Each region multiplies all branch lengths by its designed relative rate. The default WAG substitution model is independent of the scoring methods' BLOSUM62 matrix. JTT and LG are available. Method 1 is the default, and `--method 2` selects INDELible's alternative simulation algorithm.

Insertion and deletion rates are event rates relative to an average substitution rate of 1. Fixed event rates default to zero so the A–D gap-free controls retain their original design. Each event type has an independent length distribution. The default is truncated Zipf, `POW 1.7 15`, to favor shorter events while bounding the tail. Supported syntax includes `NB q r`, `POW exponent [maximum]`, `LAV exponent maximum` and `USER weight1 weight2 ...`. USER weights correspond to lengths 1, 2, and so on.

With target gap occupancies, the generator tests up to seven candidate total event rates on five independently seeded pilots per candidate. It scales insertion and deletion together, keeping the supplied ratio, or using 1:1 when both rates are zero. It uses the closest tested value and reuses it for matching regions within the family. The production regions each use separate seeds. The bounded search tests total rates between 0 and 1.5. Metadata records the target, candidate measurements, selected rates and production occupancy. Finite families, short regions and the upper search bound can leave residual differences from the target.

By default the designed relative rates are 0.06, 0.32 and 1.0. Every simulated column, including insertions, inherits its source region's rate. The existing truth TSV keeps one-based `position`, `rate`, `conservation`, `conserved` and `boundary_before`. The continuous design target remains `exp(-effective_divergence × rate)` and `conserved` is 1 for rates at or below 0.18. Partition IDs, rate classes, boundary masks and observed column gap fractions are additional fields. Tree-based divergence and observed pairwise identity are reported separately from the designed target.

Separate partitions constrain indel events to their own region. Rates are constant within each partition, and the reference `gradual` scenario uses up to twelve piecewise-constant intervals. The `redundant` scenario uses a near-tip clade. The optional `--boundary-width` marks a number of columns on either side of each join in the truth TSV for sensitivity analyses. The current benchmark includes those columns under its existing coverage rules.

## Generate one family

```powershell
python -m data_generation --scenario three-rate-blocks --size 4 --length 60 --seed 42 --insertion-rate 0.08 --deletion-rate 0.05 --boundary-width 2 --output-dir data_generation/datasets/pilot
```

The command also accepts a JSON object with underscore-separated option names. For example, save these settings as a file named `pilot.json` and set the output path and seed on the command line.

```json
{
  "scenario": "three-rate-blocks",
  "size": 4,
  "length": 60,
  "divergence": 0.7,
  "seed": 42,
  "model": "WAG",
  "insertion_rate": 0.08,
  "deletion_rate": 0.05,
  "boundary_width": 2
}
```

```powershell
python -m data_generation --config pilot.json --seed 43 --output-dir data_generation/datasets/pilot_43
```

Paths supplied in JSON are relative to the configuration file. Command-line arguments override fields in the JSON object, and command-line paths are relative to the working directory. Use `--help` for other settings, including method, size, ancestral length, model, divergence, redundancy, per-event indel rates and lengths, and target gap occupancy. The default single-family design is 50 taxa, 240 ancestral residues, a root-to-tip divergence of 0.7 and seed 1.

For arbitrary regions, pass `--regions` with a JSON array ordered from N-terminus to C-terminus. Each entry supplies an ancestral `length` and relative `rate`, plus optional `gap_target` and overrides of the global settings. For example, this configuration uses `--insertion-rate` and `--deletion-rate` for the first two regions and an independent deletion rate and distribution for the final region.

```json
[
  {"length": 20, "rate": 0.06},
  {"length": 20, "rate": 0.32},
  {"length": 20, "rate": 1.0, "deletion_rate": 0.1, "deletion_length": {"name": "USER", "parameters": [6, 3, 1]}}
]
```

Choose a new output directory for each run. It contains `alignment.fasta`, `ground_truth.tsv`, `unaligned.fasta`, `metadata.json` and a raw directory for each region. Use `--output` and `--truth` for individual file paths and `--work-dir` to set the raw directory (defaults to the FASTA name with a `.run` suffix). The raw output retains case-sensitive insertion markers, ancestral FASTA, site-rate reports, simulator logs and the exact `control.txt`. Metadata includes requested ancestry and final lengths, trees, seeds, actual insertion and deletion counts and mean event lengths, aligned gap-run length counts, and mean pairwise identity and p-distance over positions with residues present in both taxa. Use the same executable build, options and seed to reproduce the sequence and truth files.

## Generate the full supplemental suites

```powershell
python -m data_generation --suite abc --output-dir data_generation/datasets/indelible_abc
python -m data_generation --suite d --output-dir data_generation/datasets/indelible_d
python -m data_generation --suite d-gaps --output-dir data_generation/datasets/indelible_d_gaps
```

| Pattern | Ancestral design | Varied setting |
| --- | --- | --- |
| A | Uniform low rate | Root-to-tip divergence 0.35, 0.7, 1.4 |
| B | Low-rate flanks and a central high-rate block | High-rate fraction 25%, 40%, 60% |
| C | A low-rate motif in a high-rate background | Motif length 3, 5, 9 residues |
| D | Low, medium, high, medium, low | High-rate fraction 25%, 40%, 60% |

The default is 50 sequences, 240 ancestral sites and four replicates per setting. The manifest's `length` is the actual alignment width. `ancestral_length` records the original length separately. The `d-gaps` suite defaults to gap occupancy targets of 3%, 12% and 25% for low, medium and high regions. The production seed is the manifest family seed minus 12, which pairs the default D-with-gaps run with the gap-free D suite. Use the same offset when overriding both suites' starting seeds. Evolutionary indels change the subsequent simulation and its alignment coordinates, so comparisons are at the family or setting level. Pattern A contains a single designed class, which makes binary AUROC and AUPRC undefined. The manifest and benchmark family table give the reason. The suite's `validation.json` reports the observed p-distance by rate and whether it increases across rates at the suite level.

Each suite writes a manifest at its root. Its FASTA, truth and metadata paths point to the same per-family layout used by single-family runs in the suite's runs directory. See the [benchmark commands](../benchmarks/README.md) for evaluation.

## Reproduce the recorded INDELible benchmark inputs

The [SCAW report](../rvne/results/scaw_revision/README.md) identifies 12 INDELible families, seeds 81101–81106 and 81201–81206. Its companion `fresh_inputs.zip` contains the FASTA, per-column truth, metadata and native controls for each seed. To regenerate a family, use the recorded tree design, options and executable build. For example, the 81106 indel family uses these settings.

```powershell
python -m data_generation --scenario three-rate-blocks --size 20 --length 120 --seed 81106 --insertion-rate 0.03 --deletion-rate 0.03 --output-dir data_generation/datasets/indelible_81106
```

Fletcher, W. and Yang, Z. (2009). [INDELible: A Flexible Simulator of Biological Sequence Evolution](https://doi.org/10.1093/molbev/msp098). *Molecular Biology and Evolution*, 26(8), 1879–1888.

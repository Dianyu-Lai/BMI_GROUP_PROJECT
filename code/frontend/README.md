# Protein Conservation Scorer

Run the demo from the repository root with Python 3.10 or newer.

```powershell
python -m pip install -r code/rvne/requirements.txt
python code/frontend/server.py
```

Open [http://127.0.0.1:8765/](http://127.0.0.1:8765/). Use `--port 8766` if the default port is busy. Scoring, bundled examples and recorded benchmark tables work offline.

## Generate and score your own test alignment

Open **Generate MSA** from **Analyze**. The server runs INDELible through `data_generation.rvne_simulation`.

The ready-to-run defaults match the generator: scenario `blocks`, 50 sequences, 240 ancestral sites, seed 1, divergence 0.7, WAG and zero insertion/deletion rates. All eleven built-in scenarios are available under their existing names. **Advanced settings** exposes WAG/LG/JTT, insertion/deletion rates, variable fraction for block scenarios, motif length for motif scenarios, redundancy for `redundant`, and optional target gap occupancy for `indels`. Blank optional values keep generator behavior: redundancy 0.7 for `redundant`; `indels` calibrates to 20% gap occupancy when both rates are zero. Supplying a gap target calibrates the indel rates to that target. Calibration runs extra pilot simulations and takes longer. The `divergent` scenario multiplies the entered divergence by 2.5; the result summary reports effective divergence.

Both browser and server validate parameters. This local demo allows 2–500 sequences, 12–2,000 ancestral sites, seeds 0–2,147,483,647, divergence 0–3, indel rates 0–1, and optional gap targets 0–0.5. Variable fractions must be strictly between 0 and 1; motif length must not exceed one quarter of ancestor length. Generation also checks the final alignment against the existing 2,000-column/1 MiB scoring limits; if insertions exceed them, reduce ancestor length, divergence or insertion rate. Simulator timeout remains the existing 60 seconds per process. Other engine options retain the generator's defaults.

After **Generate MSA**, the page shows sequence count, ancestor length, actual aligned columns and a short FASTA preview. Download **Aligned FASTA** and **Site truth TSV**, or choose **Use this MSA** to load the full FASTA in Analyze's **Aligned FASTA** mode. Scoring starts only when you click **Score alignment**. The truth export is unchanged: one row per final alignment column, with 1-based `position` including inserted columns, matching `scores.tsv`; it also retains region, rate, boundary and observed gap fields. Ancestor coordinates are not alignment coordinates.

INDELible discovery is unchanged: `INDELIBLE_BIN`, the checkout's `.tools/indelible/indelible.exe` on Windows (or `indelible` on other systems), then `PATH`. A configured path overrides automatic lookup. The page reports missing-executable/configuration errors and simulation failures. No new dependency is needed. Generation and scoring share the existing single active job and in-memory result cache (12 results, 30-minute expiry); restarting the server clears downloads. Temporary raw simulator files are removed after each request. Download the FASTA/truth before they expire; an already loaded FASTA can still be scored.

For a command-line example, use the same [INDELible setup](../data_generation/README.md#install-and-configure) and run the following from the repository root.

```powershell
python code/data_generation --check
python code/data_generation --scenario three-rate-blocks --size 20 --length 120 --seed 42 --insertion-rate 0.05 --deletion-rate 0.05 --output code/data_generation/datasets/assessor_42.fasta --truth code/data_generation/datasets/assessor_42.truth.tsv
```

Start the server, select **Aligned FASTA**, and upload `code/data_generation/datasets/assessor_42.fasta`. Choose **SCAW-RVNE**, run the analysis, then download `scores.tsv`. Compare its `position` column with `assessor_42.truth.tsv` to see how the scores track the designed regions. For a second run, change the seed and both output filenames. Run `python code/data_generation --help` to see all generation options.

## Analyze

- Choose **Aligned FASTA** for sequences that are already aligned. **Example dataset** offers a small practice FASTA and, when the SCAW results archive is available, five INDELible alignments with site truth. Patterns A and B illustrate uniform conservation and a central variable block. Pattern C contains a short conserved motif. The two Pattern D examples illustrate three mutation-rate tiers, with evolutionary indels in seed 81106.
- Choose **Unaligned FASTA** for sequences to align before scoring. The bundled UniProt globin example uses the existing [MAFFT and ClustalW integration](../rvne/docs/real_data_alignment.md). Install either tool on `PATH`, or point `BMI3_MAFFT` or `BMI3_CLUSTALW` to its executable before starting the server.
- Under **Advanced settings**, choose **Custom** to edit RVNE's `r`, `s1` and `s2`. SCAW-RVNE v2 exposes only **Segmentation penalty multiplier** (default 1). The actual penalty is derived from this alignment; a fixed radius of 3 is clipped at segment boundaries. The former lambda, radii, support and length fields have been removed from both UI and scoring API. Custom fields start with the previously selected preset's values, and the effective settings are available in the result download.
- Use the overview, pan buttons and zoom control to set the visible column range. Selecting a column in the profile, heatmap or alignment updates the other views and position details. The plot, heatmap and alignment support the arrow keys and Home/End. **Download** exports the alignment, site scores, settings and SCAW segments when available.
- The alignment and residue counts share a palette based on [MView P1](https://desmid.github.io/mview/manual/manual.html#colourmaps), inspired by the [SMART family alignment viewer](https://smart.embl.de/smart/show_info.pl?FORMAT=MVIEW&NAME=SH2&SHOW_ALIGNMENT=1). Hover over a legend entry to see its residue codes.

The server accepts a 1 MiB FASTA with up to 500 sequences and 2,000 aligned columns. It listens on loopback and retains up to 12 completed results in memory for 30 minutes.

## Generator checks

Run `python -m unittest discover -s code/frontend/tests -v` from the repository root. These tests check request validation, scenario defaults, missing INDELible, simulator errors, output size limits, job contention/expiry, both downloads, and real default/modified/motif generation followed by RVNE and SCAW-RVNE scoring. Real generation checks are skipped only when INDELible is unavailable. Existing generator regression tests can be run with `$env:PYTHONPATH = 'code;code/rvne/src'; python -m unittest discover -s code/rvne/tests -p test_indelible_simulation.py -v` on Windows.

## Benchmarks

The Benchmarks page uses the 12 INDELible alignments for the synthetic method comparison and ablations. Ablation results can be filtered by scenario, method and metric. The shown metrics can be computed directly from scores, coverage or detected boundaries without a fitted decision threshold. The P02144 globin profiles illustrate sensitivity to aligner choice, followed by context-only timings. Input fingerprints, baseline reproduction, study parameters and the archived cohort results are in the [SCAW study report](../rvne/results/scaw_revision/README.md).

The repository includes paired MAFFT and ClustalW scores and their settings in `code/rvne/results/globins/`. The page can display these results offline. To reproduce them, follow the [globin workflow](../rvne/docs/real_data_alignment.md) with MAFFT L-INS-i and `--settings results/scaw_revision/parameters.json`, writing to a new directory.

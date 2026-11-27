# Protein Conservation Scorer

Run the local Demo from the repository root with Python 3.10+.

```powershell
python -m pip install -r code/rvne/requirements.txt
python code/frontend/server.py
```

Open [http://127.0.0.1:8765/](http://127.0.0.1:8765/). Use `--port 8766` if the default port is busy. The browser UI uses native JavaScript modules and local CSS. The server requires NumPy through the existing RVNE requirements. Once installed, the UI, synthetic examples and benchmark tables run offline.

## Analyze

- Choose **Aligned FASTA** for sequences that are already aligned. **Example dataset** offers a small practice FASTA and five synthetic alignments reproduced from the historical benchmark seeds. Patterns A and B illustrate uniform conservation and a central variable block. Pattern C contains a short conserved motif. The two Pattern D examples share three mutation-rate tiers and use the same substitution seed, with rate-dependent gaps added to one of them.
- Choose **Unaligned FASTA** for sequences to align before scoring. The bundled UniProt globin example uses the existing [MAFFT and ClustalW integration](../rvne/docs/real_data_alignment.md). Install either tool on `PATH`, or point `BMI3_MAFFT` or `BMI3_CLUSTALW` to its executable before starting the server.
- Under **Advanced settings**, choose **Custom** to edit RVNE's `r`, `s1` and `s2`. SCAW-RVNE also exposes segmentation penalty, window smoothing, window radii, minimum support and segment lengths. Custom fields start with the previously selected preset's values, and the effective settings are available in the result download.
- Use the overview, pan buttons and zoom control to set the visible column range. Selecting a column in the profile, heatmap or alignment updates the other views and position details. The plot, heatmap and alignment support the arrow keys and Home/End. **Download** exports the alignment, site scores, settings and SCAW segments when available.

The server accepts a 1 MiB FASTA with up to 500 sequences and 2,000 aligned columns. It listens on loopback and retains up to 12 completed results in memory for 30 minutes.

## Benchmarks

The browser reads the tracked supplemental and reference result TSVs through the local server. Pattern C shows motif recall by width. Pattern D reports medium-rate detection, binary variable-site recall and Spearman. For the optional alignment sensitivity comparison, follow the [globin workflow](../rvne/docs/real_data_alignment.md) to generate both aligners' outputs in `code/rvne/results/globins/`, then reload **Benchmarks**.

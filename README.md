# BMI_GROUP_PROJECT

BMI3 Algorithms ICA Group 1 scores residue conservation in globular protein multiple sequence alignments. The main implementation reproduces the original relative von Neumann entropy (RVNE) score and evaluates segmentation-constrained adaptive windows (SCAW-RVNE) against published scoring methods on controlled synthetic alignments.

- [RVNE and SCAW-RVNE](code/rvne/README.md) provides installation, scoring and benchmark commands.
- [Benchmark results](code/rvne/docs/benchmark_report.md) documents held-out comparisons, ablations, timing and limitations.
- [Software documentation](code/rvne/docs/software_documentation.tex) gives the input and output specification and a worked scoring example.
- [Earlier baseline experiments](code/baselines/README.md) covers entropy and substitution-matrix baselines and the exploratory AAindex work.

"""Draw reproducible SVGs from the held-out summaries and a fixed motif example."""

import argparse
import csv
from html import escape
from pathlib import Path

import numpy as np

from data_generation.rvne_simulation import simulate
from msa_conservation.scoring.rvne import site_profiles
from msa_conservation.scoring.rvne_optimized import score_optimized
from msa_conservation.scoring.scaw import contextual_scores


INK = "#19263a"
MUTED = "#53657b"
GRID = "#dfe6f0"
BLUE = "#166b9a"
TEAL = "#087a66"
COPPER = "#a95825"


def label(x: float, y: float, text: str, *, size: int = 14, color: str = INK,
          anchor: str = "start", bold: bool = False) -> str:
    weight = 700 if bold else 400
    return (f'<text x="{x:g}" y="{y:g}" fill="{color}" text-anchor="{anchor}" '
            f'font-family="Segoe UI, Arial, sans-serif" font-size="{size}" '
            f'font-weight="{weight}">{escape(text)}</text>')


def line(x1, y1, x2, y2, color=GRID, width=1) -> str:
    return (f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{color}" stroke-width="{width}"/>')


def svg(path: Path, width: int, height: int, body: list[str]) -> None:
    """Write an SVG with a fixed background and scale-independent text."""
    start = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">'
             f'<rect width="{width}" height="{height}" fill="white"/>')
    path.write_text(start + "\n" + "\n".join(body) + "\n</svg>\n", encoding="utf-8")


def ablations(folder: Path) -> None:
    methods = ((0, "Original RVNE"), (1, "Cached RVNE"), (2, "Uniform window"),
               (3, "Gaussian window"), (5, "Segmentation only"), (6, "SCAW-RVNE"),
               (7, "Retuned RVNE"), (8, "Retuned SCAW-RVNE"))
    with (folder / "test_summary.tsv").open(newline="", encoding="utf-8") as handle:
        stats = {(row["method"], row["metric"]): float(row["mean"])
                 for row in csv.DictReader(handle, delimiter="\t") if row["scenario"] == "all"}
    body = [
        label(26, 35, "Conservation scoring on held-out families", size=24, bold=True),
        label(26, 58, "Mean across eligible held-out families", color=MUTED),
    ]
    for j, (metric, title, lower) in enumerate((
        ("spearman", "Spearman correlation", .8), ("auprc", "Area under precision–recall", .9)
    )):
        left = 230 + 301 * j
        width = 230
        body.append(label(left, 86, title, size=15, bold=True))
        for tick in np.linspace(lower, 1, 5):
            x = left + width * (tick - lower) / (1 - lower)
            body.extend((line(x, 101, x, 393), label(x, 413, f"{tick:.2f}", anchor="middle")))
        for row, (i, name) in enumerate(methods):
            y = 116 + 32 * row
            if j == 0:
                body.append(label(26, y + 5, f"E{i}  {name}", bold=i in (5, 8)))
            value = stats[f"E{i}", metric]
            x = left + width * (value - lower) / (1 - lower)
            color = TEAL if i == 8 else (COPPER if i == 5 else BLUE)
            body.extend((line(left, y, x, y, color, 3),
                         f'<circle cx="{x:.2f}" cy="{y}" r="4" fill="{color}"/>',
                         label(left + width + 12, y + 5, f"{value:.3f}")))
    svg(folder / "ablation_scores.svg", 840, 436, body)


def motif_example(folder: Path, generate=simulate) -> None:
    sample = generate(scenario="motifs", seed=30103)
    sites = site_profiles(sample.alignment)
    raw = score_optimized(sites)
    profiles = (
        ("E0  Original RVNE", raw * sites.coverage, MUTED),
        ("E2  Uniform window", contextual_scores(raw, sites.coverage, "uniform").scores, COPPER),
        ("E6  SCAW-RVNE", contextual_scores(raw, sites.coverage, "scaw").scores, TEAL),
    )
    start, stop = 101, 146
    columns = np.arange(start, stop)
    left, top, width, height = 72, 104, 724, 214
    x = lambda position: left + (position - start) / (stop - start - 1) * width
    shown = np.stack([scores[start-1:stop-1] for _, scores, _ in profiles])
    low = np.floor(shown.min() * 10) / 10
    high = np.ceil(shown.max() * 10) / 10
    y = lambda value: top + (high - value) / (high - low) * height
    body = [
        label(26, 35, "A five-column motif in a variable region", size=24, bold=True),
        label(26, 58, "Held-out synthetic family 30103", color=MUTED),
        label(26, 91, "Score", bold=True),
    ]
    truth = np.flatnonzero(sample.conserved[start-1:stop-1]) + start
    motif_left, motif_right = x(truth.min() - .5), x(truth.max() + .5)
    body.append(f'<rect x="{motif_left:.2f}" y="{top}" width="{motif_right-motif_left:.2f}" '
                f'height="{height}" fill="#e3f4eb"/>')
    for tick in np.linspace(low, high, 5):
        body.extend((line(left, y(tick), left + width, y(tick)),
                     label(left - 10, y(tick) + 5, f"{tick:.2f}", anchor="end")))
    for tick in range(105, 146, 5):
        body.append(label(x(tick), top + height + 22, str(tick), anchor="middle"))
    body.append(label(left + width / 2, 364, "Alignment position", anchor="middle"))
    for i, ((name, _, color), values) in enumerate(zip(profiles, shown)):
        points = " ".join(f"{x(column):.2f},{y(score):.2f}" for column, score in zip(columns, values))
        body.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.3" '
                    'stroke-linecap="round" stroke-linejoin="round"/>')
        sample_x = 73 + 239 * i
        body.extend((line(sample_x, 391, sample_x + 21, 391, color, 3),
                     label(sample_x + 29, 396, name)))
    body.append(label((motif_left + motif_right) / 2, 92, "True motif", color=TEAL, anchor="middle"))
    svg(folder / "short_motif.svg", 840, 418, body)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    try:
        motif_example(args.results, simulate)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    ablations(args.results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

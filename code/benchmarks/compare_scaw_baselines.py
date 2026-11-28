"""Score the archived fresh SCAW benchmark inputs with column-only baselines."""

import argparse
import csv
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from baselines.src.msa_conservation.msa import Alignment as BaselineAlignment
from baselines.src.msa_conservation.scoring.baseline import score_alignment
from baselines.src.msa_conservation.scoring.substitution import blosum62_integer_scores
from msa_conservation.msa import henikoff_weights, read_fasta
from msa_conservation.scoring.rvne import RVNEParameters, score_reference, site_profiles
from .run_rvne_benchmark import summarize, write_tsv
from .run_supplemental_benchmark import read_tsv
from .rvne_evaluation import comparator_profiles, site_metrics


COHORT = "fresh_indelible"
METRICS = ("auprc", "spearman", "auroc", "scored_fraction")
COMPARATOR_IDS = {
    "Shannon": "shannon",
    "KL": "kl",
    "JSD": "jsd",
    "BLOSUM62": "blosum62_cost",
    "ClustalW consensus": "clustalw_consensus",
    "AL2CO variance": "al2co_variance",
}


def standalone_profiles(alignment, sites):
    """Use the baseline package's distinct-pair and adjusted-entropy scores."""
    baseline_alignment = BaselineAlignment(alignment.ids, alignment.sequences)
    columns = tuple(score_alignment(baseline_alignment, weights=henikoff_weights(alignment),
                                    reference_baselines=True))
    matrix = blosum62_integer_scores()
    minimum, maximum = min(map(min, matrix)), max(map(max, matrix))
    negative_sop = np.array([column.blosum62_sop_variability for column in columns], dtype=float)
    entropy = np.array([column.adjusted_entropy_variability for column in columns], dtype=float)
    # Orient both as nonnegative conservation scores before applying the shared gap penalty.
    scores = {
        "blosum62_sop": (-negative_sop - minimum) / (maximum - minimum),
        "adjusted_entropy": 1 / (1 + entropy),
    }
    return {method: np.where(sites.valid, profile * sites.coverage, np.nan)
            for method, profile in scores.items()}


def compare(results):
    settings = json.loads((results / "parameters.json").read_text(encoding="utf-8"))
    parameters = RVNEParameters(**settings["rvne"])
    recorded = {row["seed"]: row for row in read_tsv(results / "families.tsv")
                if row["cohort"] == COHORT and row["method"] == "v2"}
    inputs = {row["seed"]: row for row in read_tsv(results / "inputs.tsv") if row["cohort"] == COHORT}
    if not recorded or inputs.keys() != recorded.keys():
        raise ValueError("Fresh family results and recorded inputs must contain the same seeds.")
    rows = []
    with ZipFile(results / "fresh_inputs.zip") as archive:
        for seed, fingerprint in inputs.items():
            alignment = read_fasta(archive.read(f"{seed}/alignment.fasta").decode("ascii").splitlines())
            truth = list(csv.DictReader(StringIO(archive.read(f"{seed}/truth.tsv").decode("ascii")), delimiter="\t"))
            rates = np.array([float(site["rate"]) for site in truth])
            if ([int(site["position"]) for site in truth] != list(range(1, alignment.length + 1))
                    or sha256("\n".join(alignment.sequences).encode()).hexdigest() != fingerprint["alignment_sha256"]
                    or sha256(rates.tobytes()).hexdigest() != fingerprint["rates_sha256"]):
                raise ValueError(f"Archived alignment or truth differs from the recorded input for seed {seed}.")

            sites = site_profiles(alignment)
            design = np.array([float(site["conservation"]) for site in truth])
            labels = np.array([site["conserved"] == "1" for site in truth])
            scores = {
                "rvne": score_reference(sites, parameters) * sites.coverage,
                **{COMPARATOR_IDS[method]: profile
                   for method, profile in comparator_profiles(sites, alignment).items()},
                **standalone_profiles(alignment, sites),
            }
            for method, profile in scores.items():
                # Only threshold-independent metrics are exported.
                values = site_metrics(profile, design, labels, 0.)
                if (not np.array_equal(np.isfinite(profile), sites.valid)
                        or values["scored_fraction"] != float(recorded[seed]["scored_fraction"])):
                    raise ValueError(f"Scored columns for {method} differ from SCAW for seed {seed}.")
                rows.append({"method": method, "scenario": recorded[seed]["scenario"], **values})
    return [{"cohort": COHORT, **row} for row in summarize(rows, METRICS)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    write_tsv(args.results / "baseline_summary.tsv", compare(args.results))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gap_treatment.py  --  Gap handling for column-wise conservation scoring of
protein multiple sequence alignments (MSAs).

BMI3 mini-project, Project 1 (conservation score for protein MSAs).
DRAFT / REFERENCE IMPLEMENTATION. Parts of this draft were written with the
help of a large language model (Claude); the team must review, test and adapt
it before submission, and the LLM use must be declared in the documentation.

----------------------------------------------------------------------------
WHAT THIS MODULE DOES
----------------------------------------------------------------------------
Input : an aligned protein FASTA (all sequences equal length, '-' or '.' = gap)
Output: one score per alignment column (higher = more conserved), plus the
        diagnostics gap_fraction and n_residues.

A column score has two parts:
  (1) a "residue part": 1 - normalised, sequence-weighted Shannon entropy of
      the 20 amino acids. (A placeholder base score! The team's real score,
      e.g. an entropy x BLOSUM / RVNE score, can replace `entropy_conservation`
      and keep all the gap logic unchanged.)
  (2) a "gap part": how gaps are allowed to change the residue part.

Four gap treatments are implemented so that they can be benchmarked/ablated:

  mode        meaning                                                 source
  ----------  ------------------------------------------------------  -------------------
  ignore      Only residues are scored, gaps are ignored.              baseline (A)
              Known flaw: 1 residue + 9 gaps looks "fully conserved".
  symbol      Gap = 21st symbol of the entropy alphabet.               baseline (B)
              Known flaw (Valdar 2002): a mostly-gapped column can
              outscore a genuinely variable column.
  penalty     score = residue_part * (1 - g)^gamma, g = weighted       Valdar 2002 (Eq. 49,
              gap fraction. gamma = 1 is the Capra & Singh (2007)      multiplicative gap
              style penalty (score x fraction of non-gap).             prong); Capra & Singh
  threshold   Columns with g > tau get NaN ("not scored"), others      Capra & Singh code
              keep the residue part (their default cutoff is 0.3).     (gap cutoff)

Orthogonal options (can be combined with any mode):
  --terminal missing   leading/trailing gaps of a sequence (typical for
                       fragments) are treated as MISSING DATA, not as gaps.
                       (Team idea, NOT from the literature: must be validated
                       by an ablation experiment.)
  --tau-extreme X      columns with g > X are never scored (NaN). Separates
                       "not scored" from "not conserved".
  --pseudocount L      shrink residue frequencies towards a uniform background
                       when a column has very few residues (small-sample guard).
  --no-weights         use uniform sequence weights instead of Henikoff weights.
  --norm valdar|ln20   entropy normalisation: ln(min(N,K)) as in Valdar 2002
                       (Eq. 51) or the fixed ln(20) of the plain Shannon score.

All default parameter values below are PLACEHOLDERS. The course feedback asks
that every parameter is justified: choose gamma / tau / tau_extreme /
pseudocount by a sensitivity analysis on simulated data, not by hand.

Complexity: O(N * L * 20) time, O(N * L) memory (N sequences, L columns).

References (to cite in the report):
  Valdar WSJ (2002) Proteins 48:227-241.                     [C_trident, gaps]
  Henikoff S, Henikoff JG (1994) J Mol Biol 243:574-578.      [position-based weights]
  Capra JA, Singh M (2007) Bioinformatics 23:1875-1882.       [gap penalty / cutoff]

Usage examples:
  python gap_treatment.py --selftest
  python gap_treatment.py -i aligned.fasta -o scores.tsv --mode penalty --gamma 1
  python gap_treatment.py -i aligned.fasta -o all_modes.tsv --compare
"""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Tuple

import numpy as np

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------
AA = "ACDEFGHIKLMNPQRSTVWY"          # the 20 standard amino acids
GAP_CHARS = "-."                      # characters treated as alignment gaps
GAP, UNKNOWN, TERMINAL = -1, -2, -3   # integer codes (residues are 0..19)
MODES = ("ignore", "symbol", "penalty", "threshold")


# --------------------------------------------------------------------------
# 1. Reading and encoding the alignment
# --------------------------------------------------------------------------
def read_fasta(path: str) -> Tuple[List[str], List[str]]:
    """Read an aligned FASTA file. Returns (ids, upper-case sequences).

    MAFFT writes lower-case letters by default, hence the .upper().
    All sequences must have the same length.
    """
    ids: List[str] = []
    chunks: List[List[str]] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()
                ids.append(header[0] if header else "")
                chunks.append([])
            else:
                if not chunks:
                    raise ValueError("FASTA file must start with a '>' header line")
                chunks[-1].append(line)
    seqs = ["".join(c).upper() for c in chunks]
    if not seqs:
        raise ValueError("no sequences found in %s" % path)
    length = len(seqs[0])
    for name, s in zip(ids, seqs):
        if len(s) != length:
            raise ValueError(
                "sequence '%s' has length %d but the first sequence has %d: "
                "this is not an alignment" % (name, len(s), length))
    return ids, seqs


def encode(seqs: List[str], terminal: str = "gap") -> np.ndarray:
    """Convert sequences to an (N, L) int8 matrix.

    0..19 residues | GAP (-1) | UNKNOWN (-2, e.g. X/B/Z/U, treated as missing)
    | TERMINAL (-3, only if terminal == 'missing': leading/trailing gaps).
    """
    if terminal not in ("gap", "missing"):
        raise ValueError("terminal must be 'gap' or 'missing'")
    n, length = len(seqs), len(seqs[0])
    try:
        raw = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8)
    except UnicodeEncodeError:
        raise ValueError("alignment contains non-ASCII characters")
    raw = raw.reshape(n, length)

    table = np.full(256, UNKNOWN, dtype=np.int8)
    for i, a in enumerate(AA):
        table[ord(a)] = i
    for g in GAP_CHARS:
        table[ord(g)] = GAP
    m = table[raw]

    if terminal == "missing":
        is_gap = m == GAP
        non_gap = ~is_gap
        has_content = non_gap.any(axis=1)
        first = non_gap.argmax(axis=1)                       # first non-gap
        last = length - 1 - non_gap[:, ::-1].argmax(axis=1)  # last non-gap
        idx = np.arange(length)[None, :]
        term = (idx < first[:, None]) | (idx > last[:, None])
        term[~has_content, :] = True      # an all-gap sequence is all "missing"
        m[term & is_gap] = TERMINAL
    return m


# --------------------------------------------------------------------------
# 2. Sequence weights (redundancy correction)
# --------------------------------------------------------------------------
def uniform_weights(m: np.ndarray) -> np.ndarray:
    n = m.shape[0]
    return np.full(n, 1.0 / n)


def henikoff_weights(m: np.ndarray) -> np.ndarray:
    """Position-based sequence weights (Henikoff & Henikoff 1994; Valdar Eq. 53).

    In every column, a residue type shared by many sequences gives each of them
    a small weight: w_i(x) = 1 / (k_x * n_x,i), k_x = number of distinct residue
    types in column x, n_x,i = how many sequences carry sequence i's residue.
    Gaps/unknowns do not contribute. Weights are summed over columns and
    normalised to sum to 1. Needs no phylogenetic tree; O(20 * N * L).
    """
    n, length = m.shape
    counts = np.zeros((20, length))
    for a in range(20):
        counts[a] = (m == a).sum(axis=0)
    k = (counts > 0).sum(axis=0)
    w = np.zeros(n)
    for a in range(20):
        with np.errstate(divide="ignore", invalid="ignore"):
            inv = np.where(counts[a] > 0, 1.0 / (k * counts[a]), 0.0)
        w += (m == a).astype(np.float64) @ inv
    total = w.sum()
    if total <= 0:
        raise ValueError("alignment contains no standard residues")
    return w / total


# --------------------------------------------------------------------------
# 3. Column statistics
# --------------------------------------------------------------------------
def column_stats(m: np.ndarray, w: np.ndarray):
    """Per-column quantities needed by every gap treatment.

    Returns
      n_res : number of sequences with a standard residue (unweighted count)
      g     : weighted gap fraction = W_gap / (W_res + W_gap)
              (UNKNOWN and TERMINAL cells are excluded from the denominator,
              i.e. treated as missing data). NaN if the column has no
              observed cell at all.
      f     : (L, 20) weighted residue weights (NOT yet normalised)
      w_res, w_gap : weighted totals of residues / gaps in the column
    """
    length = m.shape[1]
    is_res = (m >= 0).astype(np.float64)
    is_gap = (m == GAP).astype(np.float64)
    n_res = (m >= 0).sum(axis=0)
    w_res = w @ is_res
    w_gap = w @ is_gap
    obs = w_res + w_gap
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(obs > 0, w_gap / obs, np.nan)
    f = np.zeros((length, 20))
    for a in range(20):
        f[:, a] = w @ (m == a).astype(np.float64)
    return n_res, g, f, w_res, w_gap


def residue_frequencies(f: np.ndarray, w_res: np.ndarray, n_res: np.ndarray,
                        pseudocount: float = 0.0,
                        background: Optional[np.ndarray] = None) -> np.ndarray:
    """Normalise weighted counts to frequencies over the 20 residues.

    If pseudocount > 0: p = (n * p_obs + pc * q) / (n + pc), with n the number
    of residues in the column and q a background distribution (uniform by
    default). This stops a column with a single residue from looking perfectly
    conserved. Columns without residues give NaN.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        p = f / w_res[:, None]
    p[w_res <= 0, :] = np.nan
    if pseudocount > 0:
        q = np.full(20, 1.0 / 20) if background is None else np.asarray(background, float)
        q = q / q.sum()
        n = n_res[:, None].astype(float)
        p = (n * p + pseudocount * q[None, :]) / (n + pseudocount)
        p[n_res == 0, :] = np.nan
    return p


# --------------------------------------------------------------------------
# 4. Base (placeholder) conservation score: 1 - normalised entropy
# --------------------------------------------------------------------------
def entropy_conservation(p: np.ndarray, n_ref: int, n_symbols: int = 20,
                         norm: str = "valdar") -> np.ndarray:
    """Conservation = 1 - H/H_max, H = Shannon entropy (nats) of p.

    norm='valdar': H_max = ln(min(N, K))  (Valdar 2002, Eq. 51) -- the maximum
                   can only be reached if there are at least K sequences.
    norm='ln20'  : H_max = ln(K)          (plain normalised Shannon score).
    N = number of sequences, K = number of symbols (20, or 21 with gap).
    Replace this function by the team's real score (entropy x BLOSUM, RVNE...).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        plogp = np.where(p > 0, p * np.log(p), 0.0)
    h = -plogp.sum(axis=1)                      # NaN rows stay NaN
    if norm == "ln20":
        denom = math.log(n_symbols)
    elif norm == "valdar":
        denom = math.log(min(n_ref, n_symbols))
    else:
        raise ValueError("norm must be 'valdar' or 'ln20'")
    if denom <= 0:                              # single sequence: trivially conserved
        return np.where(np.isnan(h), np.nan, 1.0)
    return 1.0 - np.clip(h / denom, 0.0, 1.0)


# --------------------------------------------------------------------------
# 5. The gap treatments
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class GapConfig:
    mode: str = "penalty"                 # ignore | symbol | penalty | threshold
    gamma: float = 1.0                    # exponent of the (1-g) penalty
    tau: float = 0.3                      # cutoff for mode='threshold'
    tau_extreme: Optional[float] = None   # never score columns with g above this
    terminal: str = "gap"                 # 'gap' or 'missing'
    use_weights: bool = True              # Henikoff weights (else uniform)
    pseudocount: float = 0.0              # small-sample shrinkage
    norm: str = "valdar"                  # entropy normalisation

    def validate(self) -> None:
        if self.mode not in MODES:
            raise ValueError("mode must be one of %s" % (MODES,))
        if self.gamma < 0:
            raise ValueError("gamma must be >= 0")
        for name in ("tau", "tau_extreme"):
            v = getattr(self, name)
            if v is not None and not (0.0 <= v <= 1.0):
                raise ValueError("%s must be in [0, 1]" % name)
        if self.pseudocount < 0:
            raise ValueError("pseudocount must be >= 0")


def score_columns(seqs: List[str], cfg: GapConfig = GapConfig()) -> Dict[str, np.ndarray]:
    """Score every column of an alignment with the chosen gap treatment.

    Returns a dict of length-L arrays:
      n_residues     residues per column
      gap_fraction   weighted gap fraction g (NaN if nothing observed)
      residue_score  conservation from residues only (before any gap rule)
      final_score    final conservation (NaN = "not scored")
    """
    cfg.validate()
    m = encode(seqs, cfg.terminal)
    n_seq = m.shape[0]
    w = henikoff_weights(m) if cfg.use_weights else uniform_weights(m)

    n_res, g, f, w_res, w_gap = column_stats(m, w)
    p20 = residue_frequencies(f, w_res, n_res, cfg.pseudocount)
    residue_score = entropy_conservation(p20, n_seq, 20, cfg.norm)

    if cfg.mode == "ignore":
        final = residue_score.copy()

    elif cfg.mode == "symbol":
        obs = w_res + w_gap
        with np.errstate(divide="ignore", invalid="ignore"):
            p21 = np.concatenate([f, w_gap[:, None]], axis=1) / obs[:, None]
        p21[obs <= 0, :] = np.nan
        final = entropy_conservation(p21, n_seq, 21, cfg.norm)

    elif cfg.mode == "penalty":
        final = residue_score * np.power(1.0 - g, cfg.gamma)

    else:  # threshold
        final = np.where(g > cfg.tau, np.nan, residue_score)

    # Rules valid for every mode: all-gap columns are never scored, and
    # extremely gappy columns can be excluded ("not scored" != "not conserved").
    final = np.where(n_res == 0, np.nan, final)
    if cfg.tau_extreme is not None:
        final = np.where(g > cfg.tau_extreme, np.nan, final)

    return {"n_residues": n_res, "gap_fraction": g,
            "residue_score": residue_score, "final_score": final}


def compare_modes(seqs: List[str], cfg: GapConfig = GapConfig()) -> Dict[str, np.ndarray]:
    """Run all four treatments with the other options held fixed (for ablations)."""
    return {mode: score_columns(seqs, replace(cfg, mode=mode))["final_score"]
            for mode in MODES}


# --------------------------------------------------------------------------
# 6. Self-test on Valdar (2002) Figure-1-style columns
# --------------------------------------------------------------------------
def _single_column(chars: str, **kw) -> float:
    """Score one column given as a string such as 'LLLL------'."""
    res = score_columns(list(chars), GapConfig(**kw))
    return float(res["final_score"][0])


def selftest(verbose: bool = True) -> None:
    columns = {
        "(j) 4xL, no gaps":     "LLLL",
        "(k) 4xL + 6 gaps":     "LLLL" + "-" * 6,
        "(m) 1xL + 9 gaps":     "L" + "-" * 9,
        "(e') 5xL + 5xD":       "L" * 5 + "D" * 5,
        "(a) 10xD":             "D" * 10,
    }
    table = {name: {mode: _single_column(col, mode=mode, gamma=1.0, tau=0.3)
                    for mode in MODES} for name, col in columns.items()}

    if verbose:
        print("%-22s" % "column" + "".join("%10s" % m for m in MODES))
        for name, row in table.items():
            print("%-22s" % name + "".join(
                "%10s" % ("NaN" if math.isnan(row[m]) else "%.3f" % row[m]) for m in MODES))

    j, k, m_, e, a = (table[n] for n in columns)

    # Valdar requirement: column (j) must be more conserved than column (k).
    assert j["penalty"] > k["penalty"] > m_["penalty"]
    assert j["symbol"] > k["symbol"]
    # A one-residue, nine-gap column must not beat a half-L/half-D column.
    assert m_["penalty"] < e["penalty"]
    # Documented failures of the naive baselines (these assertions record them):
    assert math.isclose(j["ignore"], k["ignore"])        # ignore: cannot tell (j) from (k)
    assert k["symbol"] > e["symbol"]                     # symbol: gappy beats variable
    # Threshold: gappy columns are "not scored", others keep their score.
    assert math.isnan(k["threshold"]) and not math.isnan(j["threshold"])
    # An all-gap column is never scored, whatever the mode.
    for mode in MODES:
        assert math.isnan(_single_column("-" * 10, mode=mode))
    # Larger gamma -> stronger penalty.
    assert _single_column(columns["(k) 4xL + 6 gaps"], mode="penalty", gamma=3) < k["penalty"]
    # Pseudocount lowers the score of a single-residue column.
    assert (_single_column("L" + "-" * 9, mode="ignore", pseudocount=1.0)
            < _single_column("L" + "-" * 9, mode="ignore"))
    # Sequence weighting: duplicating one sequence must not change the score
    # (single column, fixed ln(20) normalisation so N does not enter).
    base = _single_column("WQR", mode="ignore", norm="ln20")
    dup = _single_column("WQ" + "R" * 21, mode="ignore", norm="ln20")
    assert math.isclose(base, dup), (base, dup)
    # ... while without weights the duplicates DO change it.
    assert not math.isclose(
        _single_column("WQR", mode="ignore", norm="ln20", use_weights=False),
        _single_column("WQ" + "R" * 21, mode="ignore", norm="ln20", use_weights=False))

    # Terminal gaps: sequence 3 starts with two gaps (a fragment).
    toy = ["AC-D", "ACGD", "--GD"]
    g_as_gap = score_columns(toy, GapConfig(terminal="gap", use_weights=False))["gap_fraction"]
    g_missing = score_columns(toy, GapConfig(terminal="missing", use_weights=False))["gap_fraction"]
    assert math.isclose(g_as_gap[0], 1 / 3) and math.isclose(g_missing[0], 0.0)
    assert math.isclose(g_as_gap[2], 1 / 3) and math.isclose(g_missing[2], 1 / 3)  # internal gap kept
    if verbose:
        print("All self-tests passed.")


# --------------------------------------------------------------------------
# 7. Command line interface
# --------------------------------------------------------------------------
def _fmt(x: float) -> str:
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else "%.5f" % x


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Gap-aware column conservation scores for a protein MSA "
                    "(aligned FASTA). Default parameter values are placeholders.")
    ap.add_argument("-i", "--input", help="aligned FASTA file")
    ap.add_argument("-o", "--output", help="output TSV (default: stdout)")
    ap.add_argument("--mode", choices=MODES, default="penalty",
                    help="gap treatment (default: penalty)")
    ap.add_argument("--gamma", type=float, default=1.0,
                    help="exponent of (1-g) in mode=penalty (default 1)")
    ap.add_argument("--tau", type=float, default=0.3,
                    help="gap-fraction cutoff in mode=threshold (default 0.3)")
    ap.add_argument("--tau-extreme", type=float, default=None,
                    help="never score columns whose gap fraction exceeds this")
    ap.add_argument("--terminal", choices=("gap", "missing"), default="gap",
                    help="treat leading/trailing gaps as gaps or as missing data")
    ap.add_argument("--no-weights", action="store_true",
                    help="use uniform sequence weights instead of Henikoff")
    ap.add_argument("--pseudocount", type=float, default=0.0,
                    help="pseudocount towards uniform background (default 0)")
    ap.add_argument("--norm", choices=("valdar", "ln20"), default="valdar",
                    help="entropy normalisation (default valdar)")
    ap.add_argument("--compare", action="store_true",
                    help="write the final score of all four modes side by side")
    ap.add_argument("--selftest", action="store_true",
                    help="run the built-in tests and exit")
    return ap


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.selftest:
        selftest()
        return 0
    if not args.input:
        print("error: -i/--input is required (or use --selftest)", file=sys.stderr)
        return 2

    _, seqs = read_fasta(args.input)
    cfg = GapConfig(mode=args.mode, gamma=args.gamma, tau=args.tau,
                    tau_extreme=args.tau_extreme, terminal=args.terminal,
                    use_weights=not args.no_weights, pseudocount=args.pseudocount,
                    norm=args.norm)

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        if args.compare:
            base = score_columns(seqs, cfg)
            allm = compare_modes(seqs, cfg)
            out.write("pos\tn_residues\tgap_fraction\t" + "\t".join(MODES) + "\n")
            for i in range(len(base["final_score"])):
                out.write("%d\t%d\t%s\t%s\n" % (
                    i + 1, base["n_residues"][i], _fmt(float(base["gap_fraction"][i])),
                    "\t".join(_fmt(float(allm[mo][i])) for mo in MODES)))
        else:
            res = score_columns(seqs, cfg)
            out.write("pos\tn_residues\tgap_fraction\tresidue_score\tfinal_score\n")
            for i in range(len(res["final_score"])):
                out.write("%d\t%d\t%s\t%s\t%s\n" % (
                    i + 1, res["n_residues"][i], _fmt(float(res["gap_fraction"][i])),
                    _fmt(float(res["residue_score"][i])), _fmt(float(res["final_score"][i]))))
    finally:
        if args.output:
            out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())

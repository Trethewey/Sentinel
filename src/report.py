"""Assemble per-sample and per-pair score tables.

The headline column driving PASS / WARN / FAIL is `top_score_homalt`. The
combined `score` and `score_het` are retained as supplementary columns for
diagnostic review. This module emits TSVs only; HTML rendering lives in
`make_html_report.py`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

from config import RESULTS_DIR, VERDICT


# ---------------------------------------------------------------------------
# Verdict logic
# ---------------------------------------------------------------------------

def _verdict_for_sample(row: pd.Series, vd: dict = VERDICT) -> str:
    """Apply the verdict rule to one per-sample row."""
    top_homalt = row.get("top_score_homalt", np.nan)
    n_inf = row.get("top_n_informative_homalt", 0)
    deflation = row.get("homalt_deflation", np.nan)
    frac3 = row.get("frac_pairs_3plus_haps", np.nan)
    tail_frac = row.get("vaf_tail_fraction", np.nan)
    umi_mult = row.get("umi_overlap_over_background", np.nan)

    def ge(x, t):
        return (x is not None) and (not pd.isna(x)) and (x >= t)

    # FAIL conditions
    if (ge(top_homalt, vd["fail_score_homalt"]) and ge(n_inf, vd["fail_n_informative_min"])) \
       or ge(frac3, vd["fail_frac_3plus_haps"]) \
       or ge(deflation, vd["fail_homalt_deflation"]):
        return "FAIL"

    # WARN conditions (half-thresholds + softer signals)
    if ge(top_homalt, vd["warn_score_homalt"]) \
       or ge(deflation, vd["warn_homalt_deflation"]) \
       or ge(frac3, vd["warn_frac_3plus_haps"]) \
       or ge(tail_frac, vd["warn_vaf_tail_fraction"]) \
       or ge(umi_mult, vd["warn_umi_overlap_mult"]):
        return "WARN"

    return "PASS"


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def assemble_per_sample(
    samples: list,
    sample_to_group: Dict[str, int],
    n_called_sites: pd.Series,
    directional: dict,
    charr: pd.DataFrame,
    tail: pd.DataFrame,
    haplotypes: Optional[pd.DataFrame] = None,
    molecules: Optional[pd.DataFrame] = None,
    umi_jaccard: Optional[pd.DataFrame] = None,
    umi_background: Optional[float] = None,
) -> pd.DataFrame:
    """Build the per-sample report dataframe.

    All inputs are optional except samples, directional, charr, tail. Missing
    pieces simply leave their columns NaN; the verdict gate handles NaN as
    "no signal from this detector".
    """
    score_homalt = directional["score_homalt"]
    n_inf_homalt = directional["n_informative_homalt"]
    background = directional["background"]

    score_no_diag = score_homalt.where(~np.eye(len(score_homalt), dtype=bool))

    rows = []
    for s in samples:
        col = score_no_diag[s]
        # Argmax over donor column (ignore NaN); we want the donor S that
        # contributes the worst score TO this recipient T.
        if col.dropna().empty:
            top_donor = None
            top_homalt = np.nan
            top_n_inf = 0
        else:
            top_donor = col.idxmax()
            top_homalt = float(col.max())
            top_n_inf = int(n_inf_homalt.loc[top_donor, s]) if top_donor is not None else 0

        row = {
            "sample_id": s,
            "identity_group": sample_to_group.get(s),
            "n_called_sites": int(n_called_sites.get(s, 0)),
            "n_homalt_sites":      charr.loc[s, "n_homalt_sites"]      if s in charr.index else np.nan,
            "homalt_vaf":          charr.loc[s, "homalt_vaf"]          if s in charr.index else np.nan,
            "homalt_deflation":    charr.loc[s, "homalt_deflation"]    if s in charr.index else np.nan,
            "vaf_tail_fraction":   tail.loc[s, "vaf_tail_fraction"]    if s in tail.index  else np.nan,
            "top_donor_sample_id": top_donor,
            "top_score_homalt":    top_homalt,
            "top_n_informative_homalt": top_n_inf,
            "background_alt_vaf":  float(background.get(s, np.nan)),
        }

        if haplotypes is not None and s in haplotypes.index:
            row["frac_pairs_3plus_haps"]   = haplotypes.loc[s, "frac_3plus"]
            row["mean_minor_hap_fraction"] = haplotypes.loc[s, "mean_minor_hap_fraction"]
            row["n_haplotype_pairs_eval"]  = haplotypes.loc[s, "n_pairs_evaluated"]
        else:
            row["frac_pairs_3plus_haps"]   = np.nan
            row["mean_minor_hap_fraction"] = np.nan
            row["n_haplotype_pairs_eval"]  = 0

        if molecules is not None and s in molecules.index:
            row["n_consensus_reads"]         = molecules.loc[s, "n_consensus_reads"]
            row["median_molecules_per_site"] = molecules.loc[s, "median_molecules_per_site"]
            row["dup_factor"]                = molecules.loc[s, "dup_factor"]
        else:
            row["n_consensus_reads"] = np.nan
            row["median_molecules_per_site"] = np.nan
            row["dup_factor"] = np.nan

        if umi_jaccard is not None and s in umi_jaccard.index:
            off_diag = umi_jaccard.loc[s].drop(s, errors="ignore")
            if not off_diag.dropna().empty:
                row["max_umi_overlap_partner"] = off_diag.idxmax()
                row["max_umi_overlap_jaccard"] = float(off_diag.max())
                if umi_background and umi_background > 0:
                    row["umi_overlap_over_background"] = float(off_diag.max()) / umi_background
                else:
                    row["umi_overlap_over_background"] = np.nan
            else:
                row["max_umi_overlap_partner"] = None
                row["max_umi_overlap_jaccard"] = np.nan
                row["umi_overlap_over_background"] = np.nan
        else:
            row["max_umi_overlap_partner"] = None
            row["max_umi_overlap_jaccard"] = np.nan
            row["umi_overlap_over_background"] = np.nan

        rows.append(row)

    df = pd.DataFrame(rows)
    df["verdict"] = df.apply(_verdict_for_sample, axis=1)
    # Stable, readable column order
    col_order = [
        "sample_id", "identity_group", "verdict",
        "top_score_homalt", "top_donor_sample_id", "top_n_informative_homalt",
        "homalt_deflation", "homalt_vaf", "n_homalt_sites",
        "frac_pairs_3plus_haps", "mean_minor_hap_fraction", "n_haplotype_pairs_eval",
        "vaf_tail_fraction",
        "max_umi_overlap_jaccard", "max_umi_overlap_partner", "umi_overlap_over_background",
        "n_consensus_reads", "median_molecules_per_site", "dup_factor",
        "n_called_sites", "background_alt_vaf",
    ]
    return df[[c for c in col_order if c in df.columns]]


def assemble_per_pair(
    directional: dict,
    sample_to_group: Dict[str, int],
    umi_jaccard: Optional[pd.DataFrame] = None,
    score_warn: Optional[float] = None,
) -> pd.DataFrame:
    """Long-form (donor, recipient) rows for pairs with non-trivial signal."""
    score_warn = score_warn if score_warn is not None else VERDICT["warn_score_homalt"]
    sh = directional["score_homalt"]
    se = directional["score_het"]
    sc = directional["score"]
    n_inf = directional["n_informative_homalt"]

    samples = list(sh.index)
    rows = []
    for s in samples:
        for t in samples:
            if s == t:
                continue
            if sample_to_group.get(s) == sample_to_group.get(t):
                continue
            v_homalt = sh.loc[s, t]
            v_het = se.loc[s, t]
            v_sc = sc.loc[s, t]
            n = int(n_inf.loc[s, t])
            umi_j = (
                float(umi_jaccard.loc[s, t])
                if umi_jaccard is not None and s in umi_jaccard.index and t in umi_jaccard.columns
                else np.nan
            )
            flagged = (not pd.isna(v_homalt) and v_homalt >= score_warn)
            if flagged or (not pd.isna(umi_j) and umi_j > 0):
                rows.append({
                    "donor": s, "recipient": t,
                    "score_homalt": v_homalt, "score_het": v_het, "score": v_sc,
                    "n_informative_homalt": n, "umi_jaccard": umi_j,
                })
    return pd.DataFrame(rows).sort_values("score_homalt", ascending=False, na_position="last")


def write_reports(per_sample: pd.DataFrame, per_pair: pd.DataFrame, out_dir: Path = RESULTS_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    per_sample.to_csv(out_dir / "per_sample_report.tsv", sep="\t", index=False, float_format="%.6f")
    per_pair.to_csv(out_dir / "per_pair_report.tsv", sep="\t", index=False, float_format="%.6f")

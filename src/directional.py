"""Directional cross-sample contamination matrix.

For each ordered pair (S, T) with S not in T's identity group:

  informative(S, T) = sites where
                       - T is called HOM_REF                       (any alt = bleed-in/error)
                       - S is HET or HOM_ALT                       (S could donate alt reads)
                       - depth in T >= min_depth_recipient
  raw_score(S, T)    = mean alt_VAF in T over informative sites
  background(T)      = mean alt_VAF in T over sites where T is HOM_REF
                       (T's own per-base error/noise baseline)
  score(S, T)        = raw_score - background

Split by S genotype to expose alpha cleanly:
  score_homalt: at S-hom_alt sites, E[alt VAF in T] = alpha + error
  score_het:    at S-het sites,     E[alt VAF in T] = alpha/2 + error
  score:        het + hom_alt mix (power, not a direct alpha estimate)

score_homalt is the headline metric since alpha (not alpha/2) is the
expected signal at donor-hom_alt / recipient-hom_ref sites, giving a
directly interpretable contamination fraction.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np
import pandas as pd

from config import GT_HET, GT_HOM_ALT, GT_HOM_REF


def directional_matrix(
    gt: pd.DataFrame, alt: pd.DataFrame, dep: pd.DataFrame,
    min_depth_recipient: int,
    sample_to_group: Optional[Dict[str, int]] = None,
) -> dict:
    samples = list(gt.columns)
    n = len(samples)

    if sample_to_group is None:
        sample_to_group = {s: i for i, s in enumerate(samples)}
    group_arr = np.array([sample_to_group[s] for s in samples])

    gt_arr = gt.to_numpy()
    alt_arr = alt.to_numpy()
    dep_arr = dep.to_numpy()
    is_het = gt_arr == GT_HET
    is_homalt = gt_arr == GT_HOM_ALT
    is_alt_bearing = is_het | is_homalt
    is_homref = gt_arr == GT_HOM_REF
    safe_dep = np.where(dep_arr > 0, dep_arr, 1)
    site_vaf = (alt_arr / safe_dep).astype(np.float32)

    t_homref_mask = is_homref & (dep_arr >= min_depth_recipient)
    bg = np.where(t_homref_mask, site_vaf, np.nan)
    background = np.nanmean(bg, axis=0)

    # Vectorised: each result cell is a sum over sites, computable as a matrix product.
    E = t_homref_mask.astype(np.float32)            # (n_sites, n_recipients)
    A_any    = is_alt_bearing.astype(np.float32)    # (n_sites, n_donors)
    A_het    = is_het.astype(np.float32)
    A_homalt = is_homalt.astype(np.float32)
    VE = site_vaf * E                                # zero where not eligible

    # Counts: n_inf[s, t] = sum_i A[i,s] * E[i,t] = (A.T @ E)[s, t]
    n_inf_any    = (A_any.T    @ E).astype(np.int32)
    n_inf_het    = (A_het.T    @ E).astype(np.int32)
    n_inf_homalt = (A_homalt.T @ E).astype(np.int32)

    # Sum of recipient VAF over informative sites: sum_vaf[s, t] = (A.T @ VE)[s, t]
    sum_vaf_any    = A_any.T    @ VE
    sum_vaf_het    = A_het.T    @ VE
    sum_vaf_homalt = A_homalt.T @ VE

    def _mean_minus_bg(sum_vaf, n_inf):
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(n_inf > 0, sum_vaf / np.maximum(n_inf, 1), np.nan)
        return mean - background[None, :]            # broadcast bg over donors

    score        = _mean_minus_bg(sum_vaf_any,    n_inf_any)
    score_het    = _mean_minus_bg(sum_vaf_het,    n_inf_het)
    score_homalt = _mean_minus_bg(sum_vaf_homalt, n_inf_homalt)

    # Supporting reads: sites where T-eligible AND S-alt-bearing AND alt count in T >= 1
    H = (t_homref_mask & (alt_arr >= 1)).astype(np.float32)   # (n_sites, n_recipients)
    n_supp = (A_any.T @ H).astype(np.int32)

    # Mask self-pairs and within-group pairs.
    same_group = group_arr[:, None] == group_arr[None, :]     # (n, n)
    np.fill_diagonal(same_group, True)
    score[same_group]        = np.nan
    score_het[same_group]    = np.nan
    score_homalt[same_group] = np.nan
    n_inf_any[same_group]    = 0
    n_inf_het[same_group]    = 0
    n_inf_homalt[same_group] = 0
    n_supp[same_group]       = 0

    return {
        "score":                pd.DataFrame(score,        index=samples, columns=samples),
        "score_het":            pd.DataFrame(score_het,    index=samples, columns=samples),
        "score_homalt":         pd.DataFrame(score_homalt, index=samples, columns=samples),
        "n_informative":        pd.DataFrame(n_inf_any,    index=samples, columns=samples),
        "n_informative_het":    pd.DataFrame(n_inf_het,    index=samples, columns=samples),
        "n_informative_homalt": pd.DataFrame(n_inf_homalt, index=samples, columns=samples),
        "n_supporting":         pd.DataFrame(n_supp,       index=samples, columns=samples),
        "background":           pd.Series(background, index=samples, name="background_alt_vaf"),
    }

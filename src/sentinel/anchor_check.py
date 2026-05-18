"""Sample-sheet anchor check.

For each sample, compute concordance between its observed genotype and the
genotype of its EXPECTED identity (from the sample sheet). A clean sample
matches its expected identity at >0.97. A swap / heavily contaminated /
mis-labeled sample matches a different sample better than the expected one.

Outputs three per-sample fields:
  concordance_to_expected     -- 0..1, how well this BAM matches its label
  identity_match              -- True/False (>= threshold)
  best_match_sample_id        -- if mismatch, which sample it actually matches
"""
from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd

from .config import GT_NOCALL


CONCORDANCE_THRESHOLD = 0.95


def concordance(g_i: np.ndarray, g_j: np.ndarray, dep_i: np.ndarray, dep_j: np.ndarray,
                min_depth: int) -> float:
    """Strict GT equality on jointly-called sites."""
    called = (g_i != GT_NOCALL) & (g_j != GT_NOCALL) & (dep_i >= min_depth) & (dep_j >= min_depth)
    if called.sum() < 100:
        return float("nan")
    return float((g_i[called] == g_j[called]).mean())


def anchor_check(
    gt: pd.DataFrame, dep: pd.DataFrame,
    expected_identity: Dict[str, str],
    min_depth: int = 20,
    threshold: float = CONCORDANCE_THRESHOLD,
) -> pd.DataFrame:
    """For each sample, compute concordance vs its expected identity.

    expected_identity: dict {sample_id -> expected_sample_id}.
    If expected == sample_id (self), concordance is trivially 1.0 and we
    treat as "no anchor available" → identity_match = True, best_match = self.
    """
    samples = list(gt.columns)
    g = gt.to_numpy(); d = dep.to_numpy()
    sample_idx = {s: i for i, s in enumerate(samples)}

    rows = []
    for s in samples:
        i = sample_idx[s]
        exp = expected_identity.get(s, s)
        # If expected == self, no real anchor - use the cohort to find best match instead
        if exp == s or exp not in sample_idx:
            # find best-matching other sample
            best = None
            best_c = -1.0
            for t in samples:
                if t == s:
                    continue
                c = concordance(g[:, i], g[:, sample_idx[t]],
                                d[:, i], d[:, sample_idx[t]], min_depth)
                if not np.isnan(c) and c > best_c:
                    best_c = c
                    best = t
            rows.append({
                "sample_id": s,
                "expected_identity": exp,
                "concordance_to_expected": np.nan,  # no anchor
                "identity_match": True,  # nothing to disprove
                "best_match_sample_id": best,
                "best_match_concordance": best_c if best_c >= 0 else np.nan,
            })
            continue
        j = sample_idx[exp]
        c_exp = concordance(g[:, i], g[:, j], d[:, i], d[:, j], min_depth)
        # Best match among everyone else
        best = None
        best_c = -1.0
        for t in samples:
            if t == s:
                continue
            ct = concordance(g[:, i], g[:, sample_idx[t]],
                             d[:, i], d[:, sample_idx[t]], min_depth)
            if not np.isnan(ct) and ct > best_c:
                best_c = ct
                best = t
        match = (not np.isnan(c_exp)) and c_exp >= threshold
        rows.append({
            "sample_id": s,
            "expected_identity": exp,
            "concordance_to_expected": c_exp,
            "identity_match": match,
            "best_match_sample_id": best,
            "best_match_concordance": best_c if best_c >= 0 else np.nan,
        })
    return pd.DataFrame(rows).set_index("sample_id")

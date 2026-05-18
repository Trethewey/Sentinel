"""Per-chromosome breakdown of the directional score (P2b, optional).

Only runs over the flagged pairs to avoid the n^2 x n_chrom blow-up. For a
real contamination event the score is roughly uniform across chromosomes;
score concentrated on one or two chromosomes is the fingerprint of a
subclonal somatic artefact, not contamination.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import GT_HOM_ALT, GT_HOM_REF


def per_chrom_scores(
    gt: pd.DataFrame,
    alt: pd.DataFrame,
    dep: pd.DataFrame,
    site_chroms: pd.Series,
    flagged_pairs: Iterable[Tuple[str, str]],
    min_depth_recipient: int,
) -> pd.DataFrame:
    """For each (donor, recipient) in flagged_pairs, score per chromosome.

    Restricts the directional matrix algorithm to S=hom_alt sites only
    (matching the score_homalt headline). Returns long-form DataFrame with:
        donor, recipient, chrom, score_homalt, n_informative_homalt
    """
    samples = list(gt.columns)
    g = gt.to_numpy(); a = alt.to_numpy(); d = dep.to_numpy()
    chroms = site_chroms.to_numpy()
    safe_dep = np.where(d > 0, d, 1)
    vaf = a / safe_dep
    is_homalt = g == GT_HOM_ALT
    is_homref = g == GT_HOM_REF
    sample_idx = {s: i for i, s in enumerate(samples)}
    unique_chroms = pd.unique(chroms)

    rows: List[dict] = []
    for donor, recipient in flagged_pairs:
        si = sample_idx[donor]; ti = sample_idx[recipient]
        t_homref = is_homref[:, ti] & (d[:, ti] >= min_depth_recipient)
        # T's background per chromosome
        for c in unique_chroms:
            cmask = chroms == c
            bg_mask = cmask & t_homref
            if bg_mask.sum() == 0:
                background = 0.0
            else:
                background = float(vaf[bg_mask, ti].mean())
            inf_mask = cmask & t_homref & is_homalt[:, si]
            n_inf = int(inf_mask.sum())
            if n_inf == 0:
                score = np.nan
            else:
                score = float(vaf[inf_mask, ti].mean()) - background
            rows.append({
                "donor": donor, "recipient": recipient, "chrom": c,
                "score_homalt": score, "n_informative_homalt": n_inf,
            })
    return pd.DataFrame(rows)

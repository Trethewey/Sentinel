"""CHARR-style hom-alt VAF deflation per sample (P1a).

For a pure diploid sample, the mean alt VAF at its OWN called hom-alt sites
should be ~1.0 (with a small dip from sequencing error). Any source of
non-self reads pulls that mean down. Concretely, if a fraction alpha of
reads come from a donor:

  - For sites where donor is hom-ref, contributed VAF -> 0   (deflation = alpha)
  - For sites where donor is het,     contributed VAF -> 0.5 (deflation = alpha/2)
  - For sites where donor is hom-alt, contributed VAF -> 1.0 (no deflation)

Averaged over a generic donor with ~half het / half hom-ref at this sample's
hom-alt sites, the expected per-site mean alt VAF is approximately
  1 - alpha * (P(donor_homref) + 0.5 * P(donor_het))
i.e. deflation grows linearly with alpha. This is orthogonal to the
cross-sample directional matrix and detects off-flowcell contamination
(reagents, cell-line stock, environmental) that the matrix cannot see.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import GT_HOM_ALT


def homalt_vaf_deflation(
    gt: pd.DataFrame, alt: pd.DataFrame, dep: pd.DataFrame,
    min_depth_recipient: int,
) -> pd.DataFrame:
    """Return per-sample homalt VAF deflation summary.

    Columns:
      homalt_vaf         mean alt-VAF at sample's own hom-alt sites
      homalt_deflation   1 - homalt_vaf  (pure ~ small; contaminated grows with alpha)
      n_homalt_sites     number of hom-alt sites that passed the depth filter
    """
    samples = list(gt.columns)
    g = gt.to_numpy(); a = alt.to_numpy(); d = dep.to_numpy()

    safe_dep = np.where(d > 0, d, 1)
    vaf = a / safe_dep

    rows = []
    for j, s in enumerate(samples):
        mask = (g[:, j] == GT_HOM_ALT) & (d[:, j] >= min_depth_recipient)
        n_sites = int(mask.sum())
        if n_sites == 0:
            rows.append((s, np.nan, np.nan, 0))
            continue
        mv = float(vaf[mask, j].mean())
        rows.append((s, mv, 1.0 - mv, n_sites))
    return pd.DataFrame(
        rows,
        columns=["sample_id", "homalt_vaf", "homalt_deflation", "n_homalt_sites"],
    ).set_index("sample_id")

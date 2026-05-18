"""Per-sample VAF-in-tails anomaly score.

Single-sample, cohort-independent QC: fraction of called sites whose alt VAF
falls into the "shoulder" bands [0.05, 0.15] or [0.85, 0.95] - i.e. neither
clean hom-ref/hom-alt nor clean het. Elevated values indicate contamination,
mapping issues, or library prep problems. Does not replace anything; gives
an at-a-glance per-sample dashboard metric.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
import pandas as pd

from .config import TAIL_HI, TAIL_LO


def vaf_tail_fraction(
    gt: pd.DataFrame,
    alt: pd.DataFrame,
    dep: pd.DataFrame,
    min_depth: int,
    lo: Tuple[float, float] = TAIL_LO,
    hi: Tuple[float, float] = TAIL_HI,
) -> pd.Series:
    """Per-sample fraction of well-covered sites whose VAF lies in a tail band.

    Counts ALL well-covered sites (regardless of called genotype), so the
    metric is independent of the caller's thresholds. Sample is the column
    axis of the input matrices.
    """
    samples = list(dep.columns)
    a = alt.to_numpy(); d = dep.to_numpy()
    safe_dep = np.where(d > 0, d, 1)
    vaf = a / safe_dep
    well = d >= min_depth

    out = []
    for j, s in enumerate(samples):
        mask = well[:, j]
        n = int(mask.sum())
        if n == 0:
            out.append((s, np.nan, 0))
            continue
        v = vaf[mask, j]
        in_lo = (v >= lo[0]) & (v <= lo[1])
        in_hi = (v >= hi[0]) & (v <= hi[1])
        frac = float((in_lo | in_hi).mean())
        out.append((s, frac, n))
    return pd.DataFrame(
        out, columns=["sample_id", "vaf_tail_fraction", "n_sites_tail_eval"]
    ).set_index("sample_id")

"""Wide site-by-sample matrix construction.

Site key is 'chrom:pos:ref>alt' so the alt allele is implicitly fixed for
all downstream comparisons (important for the LOH-tolerant concordance rule).
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
import pandas as pd

from .config import GT_NOCALL


def pivot_matrices(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (gt, alt_depth, depth) wide matrices indexed by site, columns by sample.

    'site' is the concatenated key 'chrom:pos:ref>alt'. Categorical encoding is
    used so the index build does not walk Python-level string concatenations
    across millions of rows.
    """
    chrom = df["chrom"].astype(str).to_numpy()
    pos = df["pos"].astype(str).to_numpy()
    ref = df["ref"].to_numpy()
    alt = df["alt"].to_numpy()
    sites = np.char.add(np.char.add(np.char.add(np.char.add(np.char.add(
        chrom, ":"), pos), ":"), ref), np.char.add(">", alt))
    site_cat = pd.Categorical(sites)

    df = df.assign(_site=site_cat)
    df = df.set_index(["_site", "sample_id"])[["gt", "alt_depth", "depth"]]
    # One unstack does the work of three pivot_tables and preserves dtypes.
    wide = df.unstack("sample_id")

    gt_mat = wide["gt"].fillna(GT_NOCALL).astype(np.int8)
    alt_mat = wide["alt_depth"].fillna(0).astype(np.int32)
    dep_mat = wide["depth"].fillna(0).astype(np.int32)
    for m in (gt_mat, alt_mat, dep_mat):
        m.index.name = "site"
    return gt_mat, alt_mat, dep_mat


def site_chroms(matrix: pd.DataFrame) -> pd.Series:
    """Extract the chromosome from each 'chrom:pos:ref>alt' index entry."""
    idx = pd.Index(matrix.index.astype(str))
    return pd.Series(idx, index=idx).str.split(":", n=1).str[0].rename("chrom")

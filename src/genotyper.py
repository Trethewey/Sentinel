"""Per-sample genotype calling from REF/ALT allele depths."""
from __future__ import annotations

import numpy as np
import pandas as pd

from config import (
    GT_HET, GT_HOM_ALT, GT_HOM_REF, GT_NOCALL,
    MAX_VAF_HOM_REF,
)


def call_genotypes(
    ad: pd.DataFrame,
    min_depth: int,
    het_lo: float,
    het_hi: float,
    hom_alt_lo: float,
    hom_ref_hi: float = MAX_VAF_HOM_REF,
) -> pd.DataFrame:
    """Append 'vaf' and 'gt' columns to a long-form AD frame.

    Parameters
    ----------
    ad : DataFrame with at least 'depth' and 'alt_depth' columns.
    min_depth : int, minimum depth to attempt a call (otherwise GT_NOCALL).
    het_lo, het_hi : VAF bounds for a heterozygous call.
    hom_alt_lo : minimum VAF for a homozygous-alt call.
    hom_ref_hi : maximum VAF for a homozygous-ref call (default 0.05).
    """
    df = ad.copy()
    vaf = df["alt_depth"] / df["depth"].replace(0, np.nan)
    gt = np.full(len(df), GT_NOCALL, dtype=np.int8)
    enough_depth = df["depth"].to_numpy() >= min_depth
    v = vaf.to_numpy()
    gt[enough_depth & (v <= hom_ref_hi)] = GT_HOM_REF
    gt[enough_depth & (v >= het_lo) & (v <= het_hi)] = GT_HET
    gt[enough_depth & (v >= hom_alt_lo)] = GT_HOM_ALT
    df["vaf"] = vaf
    df["gt"] = gt
    return df

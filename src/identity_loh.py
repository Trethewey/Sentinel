"""LOH-tolerant identity concordance + union-find clustering.

Same-patient tumour/normal pairs frequently show loss-of-heterozygosity
in panel regions; without an LOH-tolerant rule those pairs slip below the
cluster threshold and the two samples are then scored against each other
as private donors, which manufactures false positives in the directional
matrix.

Rule (per site, restricted to sites where at least one of the two samples
is alt-bearing - sites where both are hom_ref are uninformative for
identity since most panel sites are hom_ref in most people):

    exact match (g_i == g_j)             -> concordant
    het <-> hom_alt                       -> concordant (LOH gain of alt)
    anything else                         -> discordant

The het<->hom_ref case is excluded because unrelated diploid samples have
many such sites by chance, which would collapse the cohort into one group.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from config import GT_HET, GT_HOM_ALT, GT_HOM_REF, GT_NOCALL


def _pair_masks(g_arr: np.ndarray, dep_arr: np.ndarray, min_depth: int):
    """Return per-sample boolean matrices used by both pair and matrix paths.

    Returns dict of (n_sites, n_samples) bool arrays:
      dep_ok, called, het_d, homalt_d, ab_dep
    """
    dep_ok = dep_arr >= min_depth
    called = (g_arr != GT_NOCALL) & dep_ok
    het_d = (g_arr == GT_HET) & dep_ok
    homalt_d = (g_arr == GT_HOM_ALT) & dep_ok
    ab_dep = het_d | homalt_d
    return dict(dep_ok=dep_ok, called=called, het_d=het_d, homalt_d=homalt_d, ab_dep=ab_dep)


def loh_concordance_pair(
    gt_i: pd.Series, gt_j: pd.Series,
    alt_i: pd.Series, alt_j: pd.Series,
    dep_i: pd.Series, dep_j: pd.Series,
    min_depth: int,
) -> float:
    g = np.stack([gt_i.to_numpy(), gt_j.to_numpy()], axis=1)
    d = np.stack([dep_i.to_numpy(), dep_j.to_numpy()], axis=1)
    m = _pair_masks(g, d, min_depth)
    return _concordance_from_masks(m, 0, 1)


def _concordance_from_masks(m: dict, i: int, j: int) -> float:
    called_i = m["called"][:, i]; called_j = m["called"][:, j]
    ab_i = m["ab_dep"][:, i]; ab_j = m["ab_dep"][:, j]
    both_called = called_i & called_j
    if int(both_called.sum()) < 100:
        return float("nan")
    valid = both_called & (ab_i | ab_j)
    n_valid = int(valid.sum())
    if n_valid < 50:
        return float("nan")
    # concordant at valid sites = both ab AND (same genotype or het<->homalt)
    # this equals "both ab at valid sites" since at valid sites at least one is ab,
    # and concordance categories cover exactly the both-ab cases (het-het, homalt-homalt,
    # het-homalt, homalt-het).
    n_conc = int((m["ab_dep"][:, i] & m["ab_dep"][:, j]).sum())
    return n_conc / n_valid


def identity_matrix_loh(
    gt: pd.DataFrame, alt: pd.DataFrame, dep: pd.DataFrame,
    min_depth: int,
) -> pd.DataFrame:
    """Symmetric NxN LOH-tolerant concordance matrix.

    Vectorised via boolean matrix products:
      n_concordant[i, j] = sites where both i and j are alt-bearing (with depth)
      n_valid[i, j]      = sites where both called and at least one is alt-bearing
                         = (ab_i & called_j) + (called_i & ab_j) - (ab_i & ab_j)
      n_both_called[i,j] = sites where both samples are called
    """
    samples = list(gt.columns)
    n = len(samples)

    m = _pair_masks(gt.to_numpy(), dep.to_numpy(), min_depth)
    C  = m["called"].astype(np.int32)        # (n_sites, n)
    AB = m["ab_dep"].astype(np.int32)        # (n_sites, n)

    M_bc      = C.T  @ C                     # both-called counts
    M_ab_C    = AB.T @ C                     # i ab-with-dep AND j called
    M_C_ab    = C.T  @ AB                    # i called AND j ab-with-dep
    M_ab_ab   = AB.T @ AB                    # both ab-with-dep (== concordant count)

    n_valid = M_ab_C + M_C_ab - M_ab_ab
    with np.errstate(invalid="ignore", divide="ignore"):
        conc = np.where(
            (M_bc >= 100) & (n_valid >= 50),
            M_ab_ab / np.maximum(n_valid, 1),
            np.nan,
        )
    np.fill_diagonal(conc, 1.0)
    return pd.DataFrame(conc, index=samples, columns=samples)


def cluster_identity(ident: pd.DataFrame, threshold: float = 0.95) -> dict:
    """Union-find clustering on concordance >= threshold."""
    samples = list(ident.index)
    parent = {s: s for s in samples}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    arr = ident.to_numpy()
    for i, a in enumerate(samples):
        for j in range(i + 1, len(samples)):
            v = arr[i, j]
            if not np.isnan(v) and v >= threshold:
                union(a, samples[j])

    roots = {find(s) for s in samples}
    root_to_gid = {r: gid for gid, r in enumerate(sorted(roots))}
    return {s: root_to_gid[find(s)] for s in samples}

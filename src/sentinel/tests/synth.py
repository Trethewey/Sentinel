"""Synthetic allele-depth generators for the unit tests.

A single seeded generator that builds long-form AD frames at controllable
genotype frequencies, with optional alpha-mixing of a donor sample into a
recipient. Avoids any BAM I/O.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd


def random_genotypes(n: int, rng: np.random.Generator,
                     p_het: float = 0.3, p_homalt: float = 0.1) -> np.ndarray:
    r = rng.random(n)
    g = np.zeros(n, dtype=int)
    g[r < p_homalt] = 2  # hom_alt
    g[(r >= p_homalt) & (r < p_homalt + p_het)] = 1  # het
    return g  # 0=hom_ref, 1=het, 2=hom_alt


def simulate_ad(g: np.ndarray, rng: np.random.Generator,
                depth: int = 50, err: float = 0.005,
                alpha: float = 0.0, donor_g: Optional[np.ndarray] = None):
    af = np.zeros_like(g, dtype=float)
    af[g == 1] = 0.5
    af[g == 2] = 1.0
    if alpha > 0 and donor_g is not None:
        donor_af = np.zeros_like(donor_g, dtype=float)
        donor_af[donor_g == 1] = 0.5
        donor_af[donor_g == 2] = 1.0
        af = (1 - alpha) * af + alpha * donor_af
    # sequencing-error smear
    af = np.clip(af + err * (1 - 2 * af), 0, 1)
    alt_d = rng.binomial(depth, af)
    return alt_d, np.full_like(alt_d, depth)


def long_ad_frame(genotypes_by_sample: dict, depth: int = 50,
                  rng_seed: int = 0, alpha_mixes: Optional[dict] = None) -> pd.DataFrame:
    """genotypes_by_sample: dict sample_id -> 1-D int array (0/1/2).
       alpha_mixes: optional dict sample_id -> (alpha, donor_sample_id)
                    overrides the simulated VAF for that sample.
    """
    rng = np.random.default_rng(rng_seed)
    alpha_mixes = alpha_mixes or {}
    first = next(iter(genotypes_by_sample.values()))
    n_sites = len(first)
    chroms = np.full(n_sites, "1")
    pos = np.arange(1, n_sites + 1)
    ref = np.full(n_sites, "A")
    alt = np.full(n_sites, "G")
    rows = []
    for sid, g in genotypes_by_sample.items():
        if sid in alpha_mixes:
            alpha, donor = alpha_mixes[sid]
            donor_g = genotypes_by_sample[donor]
            alt_d, d = simulate_ad(g, rng, depth=depth, alpha=alpha, donor_g=donor_g)
        else:
            alt_d, d = simulate_ad(g, rng, depth=depth)
        for i in range(n_sites):
            rows.append((sid, chroms[i], int(pos[i]), ref[i], alt[i],
                         int(d[i] - alt_d[i]), int(alt_d[i]), 0, int(d[i])))
    cols = ["sample_id", "chrom", "pos", "ref", "alt",
            "ref_depth", "alt_depth", "other_depth", "depth"]
    return pd.DataFrame(rows, columns=cols)

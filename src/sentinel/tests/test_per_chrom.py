"""Per-chromosome score breakdown tests."""
import numpy as np
import pandas as pd

from sentinel.genotyper import call_genotypes
from sentinel.per_chrom import per_chrom_scores
from sentinel.pivot import pivot_matrices, site_chroms
from sentinel.tests.synth import random_genotypes, simulate_ad


def _build_multi_chrom_frame(rng, sample_genos, depth=80, alpha_mixes=None):
    """Like synth.long_ad_frame but spreads sites across two chromosomes."""
    alpha_mixes = alpha_mixes or {}
    first = next(iter(sample_genos.values()))
    n = len(first)
    half = n // 2
    chroms = np.array(["1"] * half + ["2"] * (n - half))
    pos = np.concatenate([np.arange(1, half + 1), np.arange(1, n - half + 1)])
    ref = np.full(n, "A"); alt = np.full(n, "G")
    rows = []
    for sid, g in sample_genos.items():
        if sid in alpha_mixes:
            a, donor = alpha_mixes[sid]
            ad, d = simulate_ad(g, rng, depth=depth, alpha=a, donor_g=sample_genos[donor])
        else:
            ad, d = simulate_ad(g, rng, depth=depth)
        for i in range(n):
            rows.append((sid, chroms[i], int(pos[i]), ref[i], alt[i],
                         int(d[i] - ad[i]), int(ad[i]), 0, int(d[i])))
    return pd.DataFrame(rows, columns=[
        "sample_id", "chrom", "pos", "ref", "alt",
        "ref_depth", "alt_depth", "other_depth", "depth",
    ])


def test_uniform_contamination_scores_consistent_across_chroms():
    rng = np.random.default_rng(31)
    n = 4000
    s = random_genotypes(n, rng, p_het=0.3, p_homalt=0.15)
    d = random_genotypes(n, rng, p_het=0.3, p_homalt=0.2)
    df = _build_multi_chrom_frame(rng, {"S": s, "D": d},
                                  depth=120, alpha_mixes={"S": (0.05, "D")})
    df = call_genotypes(df, 20, 0.25, 0.75, 0.90)
    gt, alt_m, dep_m = pivot_matrices(df)
    chrom_ser = site_chroms(gt)
    out = per_chrom_scores(gt, alt_m, dep_m, chrom_ser, [("D", "S")], 20)
    assert set(out["chrom"]) == {"1", "2"}
    scores = out.dropna(subset=["score_homalt"])["score_homalt"].to_numpy()
    assert len(scores) == 2
    # Both should be positive and within 0.04 of each other (uniform contamination)
    assert all(s > 0.005 for s in scores), f"per-chrom scores too low: {scores}"
    assert abs(scores[0] - scores[1]) < 0.04, f"per-chrom scores diverge: {scores}"

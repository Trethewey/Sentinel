"""LOH-tolerant identity tests.

We construct a sample S1 and a derived sample S2 that started as a clone of
S1 then suffered LOH-retains-alt at half its heterozygous sites
(het -> hom_alt). Under strict equality these two would diverge by ~25
percent; under the LOH-tolerant rule they should cluster together.

The het <-> hom_ref case (LOH that lost the alt allele) is deliberately
NOT rescued by the concordance rule, because unrelated diploid samples have
many such sites by chance and tolerating them collapses unrelated samples
into one identity group. See identity_loh.py for the rule.
"""
import numpy as np

from sentinel.genotyper import call_genotypes
from sentinel.identity_loh import cluster_identity, identity_matrix_loh
from sentinel.pivot import pivot_matrices
from sentinel.tests.synth import long_ad_frame, random_genotypes


def _build_loh_pair(rng):
    n = 2000
    s1 = random_genotypes(n, rng, p_het=0.5, p_homalt=0.1)
    s2 = s1.copy()
    het_idx = np.flatnonzero(s2 == 1)
    rng.shuffle(het_idx)
    # 50% of het sites suffer LOH that retains the alt allele -> hom_alt.
    cut = int(0.50 * len(het_idx))
    s2[het_idx[:cut]] = 2
    return {"S1": s1, "S2": s2}


def test_loh_pair_still_concordant():
    rng = np.random.default_rng(42)
    geno = _build_loh_pair(rng)
    df = long_ad_frame(geno, depth=60, rng_seed=42)
    df = call_genotypes(df, min_depth=20, het_lo=0.25, het_hi=0.75, hom_alt_lo=0.90)
    gt, alt, dep = pivot_matrices(df)
    ident = identity_matrix_loh(gt, alt, dep, min_depth=20)
    v = ident.loc["S1", "S2"]
    assert v >= 0.95, f"LOH-tolerant concordance too low: {v:.3f}"


def test_loh_pair_clusters_together():
    rng = np.random.default_rng(7)
    geno = _build_loh_pair(rng)
    df = long_ad_frame(geno, depth=60, rng_seed=7)
    df = call_genotypes(df, min_depth=20, het_lo=0.25, het_hi=0.75, hom_alt_lo=0.90)
    gt, alt, dep = pivot_matrices(df)
    ident = identity_matrix_loh(gt, alt, dep, min_depth=20)
    groups = cluster_identity(ident, threshold=0.95)
    assert groups["S1"] == groups["S2"]


def test_unrelated_samples_stay_below_clustering_threshold():
    """Unrelated samples must stay below the 0.95 cluster threshold."""
    rng = np.random.default_rng(11)
    s1 = random_genotypes(2000, rng, p_het=0.4, p_homalt=0.15)
    s_un = rng.permutation(random_genotypes(2000, rng, p_het=0.4, p_homalt=0.15))
    df = long_ad_frame({"S1": s1, "Sunrel": s_un}, depth=60, rng_seed=11)
    df = call_genotypes(df, min_depth=20, het_lo=0.25, het_hi=0.75, hom_alt_lo=0.90)
    gt, alt, dep = pivot_matrices(df)
    ident = identity_matrix_loh(gt, alt, dep, min_depth=20)
    v = ident.loc["S1", "Sunrel"]
    groups = cluster_identity(ident, threshold=0.95)
    assert v < 0.95, f"unrelated pair concordance crossed cluster threshold: {v:.3f}"
    assert groups["S1"] != groups["Sunrel"], "unrelated samples should not cluster"

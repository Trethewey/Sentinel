"""Directional matrix tests - score_homalt should track alpha cleanly."""
import numpy as np

from directional import directional_matrix
from genotyper import call_genotypes
from identity_loh import cluster_identity, identity_matrix_loh
from pivot import pivot_matrices
from tests.synth import long_ad_frame, random_genotypes


def test_homalt_score_recovers_alpha():
    rng = np.random.default_rng(2024)
    n = 3000
    s_donor = random_genotypes(n, rng, p_het=0.3, p_homalt=0.2)
    s_recip = random_genotypes(n, rng, p_het=0.3, p_homalt=0.05)
    s_other = random_genotypes(n, rng, p_het=0.3, p_homalt=0.1)
    geno = {"D": s_donor, "R": s_recip, "U": s_other}
    df = long_ad_frame(geno, depth=80, rng_seed=2024, alpha_mixes={"R": (0.05, "D")})
    df = call_genotypes(df, 20, 0.25, 0.75, 0.90)
    gt, alt, dep = pivot_matrices(df)
    ident = identity_matrix_loh(gt, alt, dep, 20)
    groups = cluster_identity(ident, threshold=0.95)
    res = directional_matrix(gt, alt, dep, 20, sample_to_group=groups)
    sh = res["score_homalt"].loc["D", "R"]
    assert 0.02 <= sh <= 0.08, f"score_homalt(D->R) out of expected band: {sh:.4f}"
    # Unrelated donor U should not score against R
    sh_u = res["score_homalt"].loc["U", "R"]
    assert abs(sh_u) < 0.015, f"score_homalt(U->R) leaking: {sh_u:.4f}"


def test_no_contamination_yields_low_scores():
    rng = np.random.default_rng(99)
    n = 2000
    geno = {sid: random_genotypes(n, rng, p_het=0.3, p_homalt=0.1)
            for sid in ("S1", "S2", "S3")}
    df = long_ad_frame(geno, depth=60, rng_seed=99)
    df = call_genotypes(df, 20, 0.25, 0.75, 0.90)
    gt, alt, dep = pivot_matrices(df)
    ident = identity_matrix_loh(gt, alt, dep, 20)
    groups = cluster_identity(ident, threshold=0.95)
    res = directional_matrix(gt, alt, dep, 20, sample_to_group=groups)
    sh = res["score_homalt"]
    off = sh.where(~np.eye(len(sh), dtype=bool)).abs().to_numpy()
    off = off[~np.isnan(off)]
    assert off.max() < 0.01, f"clean cohort headline above 0.01: {off.max():.4f}"

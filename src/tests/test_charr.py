"""CHARR hom-alt VAF deflation tests."""
import numpy as np

from charr import homalt_vaf_deflation
from genotyper import call_genotypes
from pivot import pivot_matrices
from tests.synth import long_ad_frame, random_genotypes


def test_clean_sample_low_deflation():
    rng = np.random.default_rng(1)
    g = random_genotypes(3000, rng, p_het=0.3, p_homalt=0.2)
    df = long_ad_frame({"S": g}, depth=80, rng_seed=1)
    df = call_genotypes(df, 20, 0.25, 0.75, 0.90)
    gt, alt, dep = pivot_matrices(df)
    out = homalt_vaf_deflation(gt, alt, dep, 20)
    defl = out.loc["S", "homalt_deflation"]
    assert 0 <= defl <= 0.01, f"clean deflation should be ~0: {defl:.4f}"


def test_spiked_sample_deflation_proportional_to_alpha():
    """At recipient's hom-alt sites, donor contributes alt VAF in [0,1]
    averaging ~0.5 across mixed genotypes; so deflation ~ alpha/2."""
    rng = np.random.default_rng(2)
    n = 3000
    s = random_genotypes(n, rng, p_het=0.3, p_homalt=0.25)
    d = random_genotypes(n, rng, p_het=0.3, p_homalt=0.1)
    alpha = 0.10
    df = long_ad_frame({"S": s, "D": d}, depth=120, rng_seed=2,
                       alpha_mixes={"S": (alpha, "D")})
    df = call_genotypes(df, 20, 0.25, 0.75, 0.90)
    gt, alt, dep = pivot_matrices(df)
    out = homalt_vaf_deflation(gt, alt, dep, 20)
    defl = out.loc["S", "homalt_deflation"]
    # Expected ~ alpha * (P(donor hom_ref) + 0.5*P(donor het))
    # With (0.65, 0.30, 0.05) priors that's ~ alpha*(0.65 + 0.15) = 0.08*alpha factor
    # which for alpha=0.10 gives ~0.08; loose bounds:
    assert 0.03 <= defl <= 0.15, f"spiked deflation out of band: {defl:.4f}"
    # Strictly bigger than clean sample (sanity check)
    clean_df = long_ad_frame({"S": s}, depth=120, rng_seed=99)
    clean_df = call_genotypes(clean_df, 20, 0.25, 0.75, 0.90)
    g2, a2, d2 = pivot_matrices(clean_df)
    clean_defl = homalt_vaf_deflation(g2, a2, d2, 20).loc["S", "homalt_deflation"]
    assert defl > clean_defl + 0.02, (
        f"spiked deflation should exceed clean by >0.02 (got {defl:.4f} vs {clean_defl:.4f})"
    )

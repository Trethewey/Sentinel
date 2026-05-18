"""VAF-in-tails anomaly score tests."""
import numpy as np

from genotyper import call_genotypes
from pivot import pivot_matrices
from tails import vaf_tail_fraction
from tests.synth import long_ad_frame, random_genotypes


def test_clean_sample_low_tail_fraction():
    rng = np.random.default_rng(13)
    g = random_genotypes(3000, rng, p_het=0.3, p_homalt=0.1)
    df = long_ad_frame({"S": g}, depth=80, rng_seed=13)
    df = call_genotypes(df, 20, 0.25, 0.75, 0.90)
    gt, alt, dep = pivot_matrices(df)
    out = vaf_tail_fraction(gt, alt, dep, min_depth=20)
    tf = out.loc["S", "vaf_tail_fraction"]
    assert tf < 0.03, f"clean tail fraction should be <0.03, got {tf:.4f}"


def test_spiked_sample_elevated_tail_fraction():
    rng = np.random.default_rng(14)
    n = 3000
    s = random_genotypes(n, rng, p_het=0.3, p_homalt=0.15)
    d = random_genotypes(n, rng, p_het=0.3, p_homalt=0.15)
    df = long_ad_frame({"S": s, "D": d}, depth=120, rng_seed=14,
                       alpha_mixes={"S": (0.10, "D")})
    df = call_genotypes(df, 20, 0.25, 0.75, 0.90)
    gt, alt, dep = pivot_matrices(df)
    out = vaf_tail_fraction(gt, alt, dep, min_depth=20)
    tf = out.loc["S", "vaf_tail_fraction"]
    assert tf > 0.05, f"spiked tail fraction should be >0.05, got {tf:.4f}"

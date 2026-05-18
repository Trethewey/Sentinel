"""3-haplotype counter tests (synthetic - no BAM I/O).

We feed `score_from_counters` with hand-built per-pair Counters so the test
exercises the scoring logic without needing pysam / a fixture BAM.
The pair-enumeration helper `find_site_pairs` is exercised on a tiny site DataFrame.
"""
from collections import Counter

import pandas as pd

from sentinel.read_haplotypes import find_site_pairs, score_from_counters


def test_pure_diploid_two_haplotypes():
    # one pair, two haplotypes well above noise -> no flag
    counters = [Counter({("A", "C"): 30, ("G", "T"): 25})]
    out = score_from_counters(counters, min_support=2)
    assert out["n_pairs_evaluated"] == 1
    assert out["n_pairs_with_3plus_haps"] == 0
    assert out["frac_3plus"] == 0.0


def test_three_haplotypes_flag_pair():
    counters = [Counter({("A", "C"): 30, ("G", "T"): 25, ("A", "T"): 5})]
    out = score_from_counters(counters, min_support=2)
    assert out["n_pairs_with_3plus_haps"] == 1
    assert out["frac_3plus"] == 1.0
    # minor fraction ~ 5/60
    assert abs(out["mean_minor_hap_fraction"] - 5 / 60) < 0.01


def test_singleton_haplotypes_filtered_out():
    # The third haplotype has a single supporting molecule -> dropped
    counters = [Counter({("A", "C"): 20, ("G", "T"): 18, ("A", "T"): 1})]
    out = score_from_counters(counters, min_support=2)
    assert out["n_pairs_with_3plus_haps"] == 0


def test_score_aggregates_across_pairs():
    counters = [
        Counter({("A", "C"): 20, ("G", "T"): 18}),                      # 2 haps -> not flagged
        Counter({("A", "C"): 15, ("G", "T"): 14, ("A", "T"): 4}),      # 3 haps -> flagged
        Counter({("A", "C"): 22, ("G", "T"): 17, ("A", "T"): 3, ("G", "C"): 3}),  # 4 haps
        Counter(),  # empty: skipped entirely
    ]
    out = score_from_counters(counters, min_support=2)
    assert out["n_pairs_evaluated"] == 3
    assert out["n_pairs_with_3plus_haps"] == 2
    assert abs(out["frac_3plus"] - 2 / 3) < 1e-9


def test_alpha_estimate_from_minor_fraction():
    # With alpha=0.05, minor_hap_fraction at a 3-hap locus should be ~ 0.05
    # (both donor haplotypes carry alpha mass, but only the non-overlapping
    # one shows as a third haplotype). We test the closer-to-alpha case
    # where ONE of the two donor haplotypes is novel.
    n_recip_h1 = int(0.475 * 200); n_recip_h2 = int(0.475 * 200)
    n_donor_novel = int(0.05 * 200)
    counters = [Counter({
        ("A", "C"): n_recip_h1,
        ("G", "T"): n_recip_h2,
        ("A", "T"): n_donor_novel,
    })]
    out = score_from_counters(counters, min_support=2)
    minor = out["mean_minor_hap_fraction"]
    assert 0.03 <= minor <= 0.07, f"alpha estimate from minor_frac off: {minor:.3f}"


def test_find_site_pairs_within_distance():
    sites = pd.DataFrame([
        {"chrom": "1", "pos": 100, "ref": "A", "alt": "G"},
        {"chrom": "1", "pos": 250, "ref": "C", "alt": "T"},   # 150 bp away -> pair
        {"chrom": "1", "pos": 500, "ref": "T", "alt": "A"},   # 250 bp from 250 -> pair
        {"chrom": "1", "pos": 5000, "ref": "G", "alt": "C"},  # too far
        {"chrom": "2", "pos": 100, "ref": "A", "alt": "G"},   # different chrom
    ])
    pairs = find_site_pairs(sites, max_dist=300)
    # Expected: (100,250), (250,500) on chrom 1; nothing else
    assert len(pairs) == 2
    positions = sorted([(a[1], b[1]) for a, b in pairs])
    assert positions == [(100, 250), (250, 500)]

"""UMI Jaccard + footprint sanity tests.

Tests exercise the pure-Python functions (jaccard, umi_footprint_overlap,
background_jaccard_estimate) using synthetic UMI sets - no BAM I/O.
"""
import math

import numpy as np

from umi_metrics import (
    background_jaccard_estimate, jaccard, umi_footprint_overlap,
)


def _make_umis(prefix: str, n: int):
    return {f"{prefix}{i:06d}" for i in range(n)}


def test_jaccard_disjoint_zero():
    a = _make_umis("A", 1000)
    b = _make_umis("B", 1000)
    assert jaccard(a, b) == 0.0


def test_jaccard_identical_one():
    a = _make_umis("A", 100)
    assert jaccard(a, a) == 1.0


def test_jaccard_known_overlap():
    a = set(range(100))
    b = set(range(80, 180))
    # |inter| = 20, |union| = 180
    assert abs(jaccard(a, b) - 20 / 180) < 1e-9


def test_footprint_matrix_symmetric_and_diagonal_one():
    umis = {
        "S1": set(range(0, 1000)),
        "S2": set(range(500, 1500)),
        "S3": set(range(2000, 3000)),
    }
    m = umi_footprint_overlap(umis)
    assert m.shape == (3, 3)
    assert m.loc["S1", "S1"] == 1.0
    assert m.loc["S1", "S2"] == m.loc["S2", "S1"]
    assert m.loc["S1", "S3"] == 0.0
    # |inter(S1,S2)| = 500, |union| = 1500
    assert abs(m.loc["S1", "S2"] - 500 / 1500) < 1e-9


def test_background_estimate_matches_4_to_the_minus_L():
    for L in (6, 8, 10, 12):
        assert math.isclose(background_jaccard_estimate(L), 4.0 ** (-L))


def test_spiked_overlap_well_above_background():
    """Two samples share 1% of one's UMIs (simulating post-library
    contamination) - measured overlap should exceed background by orders
    of magnitude."""
    a = _make_umis("X", 10_000)
    b = _make_umis("Y", 9_900) | set(list(a)[:100])  # 1% of a leaks into b
    overlap = jaccard(a, b)
    bg = background_jaccard_estimate(10)  # ~1e-6
    assert overlap > 100 * bg, f"spiked overlap {overlap} not >> background {bg}"

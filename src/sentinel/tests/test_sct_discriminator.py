"""Unit tests for the SCT vs contamination discriminator."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sentinel.sct_discriminator import call, SCTCall


def test_pass_when_nuclear_low():
    r = call(nuclear_score_homalt=0.01, frac_3plus_hap=0.30, mt_frac_mixed=0.30)
    assert r.verdict == "PASS"


def test_sct_call():
    r = call(nuclear_score_homalt=0.85, frac_3plus_hap=0.001, mt_frac_mixed=0.0)
    assert r.verdict == "SCT_LIKE"
    assert "stem cell transplant" in r.rationale.lower()


def test_contamination_call():
    r = call(nuclear_score_homalt=0.40, frac_3plus_hap=0.12, mt_frac_mixed=0.05)
    assert r.verdict == "CONTAMINATION"
    assert "contamination" in r.rationale.lower()


def test_indeterminate_disagree():
    # Strong nuclear, monoclonal hap, but mixed mt -> indeterminate
    r = call(nuclear_score_homalt=0.50, frac_3plus_hap=0.001, mt_frac_mixed=0.10)
    assert r.verdict == "INDETERMINATE"
    r2 = call(nuclear_score_homalt=0.50, frac_3plus_hap=0.10, mt_frac_mixed=0.0)
    assert r2.verdict == "INDETERMINATE"


def test_intermediate_zone_indeterminate():
    # Both signals in the "between thresholds" zone
    r = call(nuclear_score_homalt=0.30, frac_3plus_hap=0.03, mt_frac_mixed=0.01)
    assert r.verdict == "INDETERMINATE"

"""End-to-end smoke test for the scoring pipeline (matrix-only branch).

Builds a synthetic AD parquet, runs `pipeline.run` with --skip-bam-features,
and asserts the report files exist with the expected columns.
"""
import numpy as np
import pytest

from sentinel.pipeline import run
from sentinel.tests.synth import long_ad_frame, random_genotypes


@pytest.fixture
def synthetic_ad_parquet(tmp_path):
    rng = np.random.default_rng(0)
    n = 1500
    geno = {sid: random_genotypes(n, rng, p_het=0.3, p_homalt=0.15)
            for sid in ("A", "B", "C")}
    df = long_ad_frame(geno, depth=80, rng_seed=0,
                       alpha_mixes={"A": (0.05, "B")})
    p = tmp_path / "ad_long.parquet"
    df.to_parquet(p, index=False)
    return p, tmp_path


def test_pipeline_runs_and_writes_reports(synthetic_ad_parquet):
    ad_path, tmp = synthetic_ad_parquet
    out_dir = tmp / "results"
    res = run(
        ad_parquet=ad_path, out_dir=out_dir,
        skip_bam_features=True, verbose=False,
    )
    per_sample = res["per_sample"]
    assert (out_dir / "per_sample_report.tsv").exists()
    assert (out_dir / "per_pair_report.tsv").exists()
    expected_cols = {
        "sample_id", "identity_group", "verdict",
        "top_score_homalt", "top_donor_sample_id", "top_n_informative_homalt",
        "homalt_deflation", "homalt_vaf",
        "vaf_tail_fraction", "background_alt_vaf",
    }
    missing = expected_cols - set(per_sample.columns)
    assert not missing, f"per-sample report missing columns: {missing}"
    # Recipient A should have non-trivial headline score
    a_row = per_sample.set_index("sample_id").loc["A"]
    assert a_row["top_score_homalt"] > 0.01, (
        f"A's top headline score too low: {a_row['top_score_homalt']:.4f}"
    )

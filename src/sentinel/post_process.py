"""Post-processing stage.

Takes the per-sample scores from `pipeline.py` and the cached allele-depth
parquet, then:
  - runs an anchor check (concordance to each sample's declared identity)
  - applies sample-type-aware verdict thresholds
  - lists tied candidate donors when scores are within a tolerance band

Paths come from the PROJ environment variable. Outputs land in <PROJ>/results/.
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from .anchor_check import anchor_check, CONCORDANCE_THRESHOLD
from .config import MIN_DEPTH_CALL
from .genotyper import call_genotypes
from .io_ad import load_ad_long
from .pivot import pivot_matrices


PROJ = Path(os.environ.get("PROJ", "."))
SCORES_DIR = PROJ / "results" / "scores"
RESULTS_DIR = PROJ / "results"
MANIFEST = PROJ / "work" / "run_manifest.tsv"
AD = PROJ / "work" / "ad_long.parquet"


# Sample-type-aware verdict thresholds. Controls don't flag (they're inherently
# spike-in or known-noisy material). FFPE tolerates more noise than fresh.
VERDICT_RULES = {
    "control":         dict(score=99.0, defl=99.0,  tail=99.0,  warn_score=99.0,  warn_defl=99.0,  warn_tail=99.0),
    "reference_ffpe":  dict(score=0.03, defl=0.05,  tail=0.20,  warn_score=0.015, warn_defl=0.03,  warn_tail=0.12),
    "reference_fresh": dict(score=0.02, defl=0.025, tail=0.10,  warn_score=0.008, warn_defl=0.012, warn_tail=0.07),
    "clinical":        dict(score=0.02, defl=0.025, tail=0.10,  warn_score=0.008, warn_defl=0.012, warn_tail=0.07),
    "mix":             dict(score=0.02, defl=0.025, tail=0.10,  warn_score=0.008, warn_defl=0.012, warn_tail=0.07),
}


def verdict_for(row, sample_type):
    rules = VERDICT_RULES.get(sample_type, VERDICT_RULES["clinical"])
    sh = float(row.get("top_score_homalt", 0) or 0)
    nf = float(row.get("top_n_informative_homalt", 0) or 0)
    df = float(row.get("homalt_deflation", 0) or 0)
    tf = float(row.get("vaf_tail_fraction", 0) or 0)
    if row.get("identity_match") is False:
        return "FAIL"
    if (sh >= rules["score"] and nf >= 30) or df >= rules["defl"] or tf >= rules["tail"]:
        return "FAIL"
    if (sh >= rules["warn_score"] and nf >= 30) or df >= rules["warn_defl"] or tf >= rules["warn_tail"]:
        return "WARN"
    return "PASS"


def top_n_donors(score_row: pd.Series, n: int = 3, tie_tol: float = 0.005) -> str:
    """Top-N candidate donors as a string, including every donor within tie_tol of the leader."""
    s = score_row.dropna().sort_values(ascending=False)
    if s.empty:
        return ""
    top = s.iloc[0]
    in_tie = s[s >= top - tie_tol]
    keep = in_tie.head(max(n, len(in_tie)))
    return "; ".join(f"{idx}({val:.3f})" for idx, val in keep.items())


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading allele depths, manifest, and per-sample scores")
    ad = load_ad_long(AD)
    manifest = pd.read_csv(MANIFEST, sep="\t")
    expected = dict(zip(manifest["sample_id"], manifest.get("expected_identity", manifest["sample_id"])))
    sample_type_map = dict(zip(manifest["sample_id"], manifest.get("sample_type", ["clinical"] * len(manifest))))
    well_map = dict(zip(manifest["sample_id"], manifest.get("well", [""] * len(manifest))))
    per_sample = pd.read_csv(SCORES_DIR / "per_sample_report.tsv", sep="\t").set_index("sample_id")
    score_homalt = pd.read_csv(SCORES_DIR / "directional_contamination_homalt_only.tsv", sep="\t", index_col=0)

    print("Calling genotypes for the anchor check")
    df = call_genotypes(ad, MIN_DEPTH_CALL, 0.25, 0.75, 0.90)
    gt_m, _alt_m, dep_m = pivot_matrices(df)

    print(f"Anchor check ({CONCORDANCE_THRESHOLD} concordance threshold)")
    anchor = anchor_check(gt_m, dep_m, expected, min_depth=MIN_DEPTH_CALL, threshold=CONCORDANCE_THRESHOLD)
    print(f"  identity matches: {int(anchor['identity_match'].sum())}/{len(anchor)}")

    print("Building tied-donor lists")
    top_n = score_homalt.apply(lambda col: top_n_donors(col), axis=0)
    top_n.name = "top_donor_tied_set"

    print("Producing final per-sample report")
    final = per_sample.join(anchor, how="left")
    final["sample_type"] = final.index.map(lambda s: sample_type_map.get(s, "clinical"))
    final["well"] = final.index.map(lambda s: well_map.get(s, ""))
    final["top_donor_tied_set"] = top_n.reindex(final.index)
    final["verdict"] = final.apply(lambda r: verdict_for(r, r["sample_type"]), axis=1)

    out_tsv = RESULTS_DIR / "per_sample_report.tsv"
    out_csv = RESULTS_DIR / "per_sample_report.csv"
    final.to_csv(out_tsv, sep="\t")
    final.to_csv(out_csv)
    print(f"Wrote {out_tsv}")

    summary = final["verdict"].value_counts().to_dict()
    print(f"\nVerdict summary: {summary}")
    by_type = final.groupby("sample_type")["verdict"].value_counts().unstack(fill_value=0)
    print("\nBy sample type:")
    print(by_type.to_string())


if __name__ == "__main__":
    main()

"""Sentinel scoring pipeline.

Reads the cached allele-depth table, calls genotypes, builds the site-by-sample
matrices, and emits the matrix-based and BAM-based scores that downstream
post-processing turns into a verdict.

Stages:
    io_ad.load_ad_long
        -> genotyper.call_genotypes
            -> pivot.pivot_matrices
                -> identity_loh.identity_matrix_loh  -> cluster_identity
                -> charr.homalt_vaf_deflation
                -> tails.vaf_tail_fraction
                -> directional.directional_matrix (uses identity groups)
                    -> per_chrom.per_chrom_scores (flagged pairs only)
    BAM-based branch (parallel-friendly per sample):
        read_haplotypes.haplotype_score
        umi_metrics.molecule_counts
    -> report.assemble + write_reports

The BAM branch is wrapped in --skip-bam-features so unit tests and dry-runs
can exercise the matrix-based detectors only.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from .charr import homalt_vaf_deflation


def _mc_remote(sid, bp, scout, sample_n_sites):
    """Module-level wrapper so ProcessPoolExecutor can pickle it."""
    from umi_metrics import molecule_counts as _mc_fn
    res = _mc_fn(bp, scout, sample_n_sites=sample_n_sites)
    res["sample_id"] = sid
    return res

from .config import (
    CACHE_DIR, GT_NOCALL, HAP_PAIR_MAX_DIST, HET_CALL_HI, HET_CALL_LO,
    HOM_ALT_LO, IDENTITY_CLUSTER_THRESHOLD, MIN_DEPTH_CALL, MIN_DEPTH_RECIPIENT,
    RESULTS_DIR, VERDICT, WORK_DIR,
)
from .directional import directional_matrix
from .genotyper import call_genotypes
from .identity_loh import cluster_identity, identity_matrix_loh
from .io_ad import load_ad_long
from .per_chrom import per_chrom_scores
from .pivot import pivot_matrices, site_chroms
from .read_haplotypes import find_site_pairs, haplotype_score
from .report import assemble_per_pair, assemble_per_sample, write_reports
from .tails import vaf_tail_fraction
from .umi_metrics import molecule_counts


def _flagged_pairs_from_directional(directional: dict, sample_to_group: Dict[str, int],
                                    threshold: float) -> list:
    sh = directional["score_homalt"]
    samples = list(sh.index)
    out = []
    for s in samples:
        for t in samples:
            if s == t or sample_to_group.get(s) == sample_to_group.get(t):
                continue
            v = sh.loc[s, t]
            if not pd.isna(v) and v >= threshold:
                out.append((s, t))
    return out


def run(
    ad_parquet: Path,
    out_dir: Path = RESULTS_DIR,
    min_depth_call: int = MIN_DEPTH_CALL,
    min_depth_recipient: int = MIN_DEPTH_RECIPIENT,
    het_lo: float = HET_CALL_LO,
    het_hi: float = HET_CALL_HI,
    hom_alt_lo: float = HOM_ALT_LO,
    identity_threshold: float = IDENTITY_CLUSTER_THRESHOLD,
    scout_sites_path: Optional[Path] = None,
    bam_manifest_path: Optional[Path] = None,
    skip_bam_features: bool = False,
    umi_length: int = 10,
    verbose: bool = True,
) -> dict:
    """End-to-end scoring run. Returns the result objects for programmatic use."""

    def log(msg: str) -> None:
        if verbose:
            print(msg, flush=True)

    out_dir.mkdir(parents=True, exist_ok=True)

    log(f"Loading allele depths from {ad_parquet}")
    ad = load_ad_long(ad_parquet)
    log(f"  {len(ad):,} rows; {ad['sample_id'].nunique()} samples")

    df = call_genotypes(ad, min_depth_call, het_lo, het_hi, hom_alt_lo)
    gt_m, alt_m, dep_m = pivot_matrices(df)
    log(f"Matrix: {gt_m.shape[0]:,} sites x {gt_m.shape[1]} samples")

    n_called = (gt_m != GT_NOCALL).sum(axis=0).rename("n_called_sites")

    log("Computing LOH-tolerant identity concordance matrix")
    ident = identity_matrix_loh(gt_m, alt_m, dep_m, min_depth_call)
    sample_to_group = cluster_identity(ident, threshold=identity_threshold)
    n_groups = len(set(sample_to_group.values()))
    log(f"  {n_groups} identity groups across {len(sample_to_group)} samples")

    log("Computing per-sample CHARR-style hom-alt deflation")
    charr_df = homalt_vaf_deflation(gt_m, alt_m, dep_m, min_depth_recipient)
    log("Computing per-sample VAF-in-tails anomaly fraction")
    tail_df = vaf_tail_fraction(gt_m, alt_m, dep_m, min_depth_call)

    log("Computing directional contamination matrix")
    directional = directional_matrix(
        gt_m, alt_m, dep_m, min_depth_recipient, sample_to_group=sample_to_group,
    )

    flagged = _flagged_pairs_from_directional(
        directional, sample_to_group, threshold=VERDICT["warn_score_homalt"],
    )
    log(f"  {len(flagged)} flagged donor-recipient pairs above WARN")
    per_chrom_df = pd.DataFrame()
    if flagged:
        per_chrom_df = per_chrom_scores(
            gt_m, alt_m, dep_m, site_chroms(gt_m), flagged, min_depth_recipient,
        )

    haplotypes_df = None
    molecules_df = None
    umi_jaccard = None
    umi_background = None
    if not skip_bam_features:
        if scout_sites_path is None or bam_manifest_path is None:
            log("BAM features skipped (no scout_sites or bam_manifest provided)")
        else:
            import os as _os
            import concurrent.futures as _cf
            N_WORKERS = max(1, (_os.cpu_count() or 4) - 2)

            try:
                log("Computing read-level 3-haplotype score per sample")
                scout = pd.read_csv(scout_sites_path, sep="\t")
                pairs = find_site_pairs(scout, max_dist=HAP_PAIR_MAX_DIST)
                log(f"  {len(pairs)} within-{HAP_PAIR_MAX_DIST}bp site pairs; "
                    f"running {N_WORKERS} parallel workers")
                manifest = pd.read_csv(bam_manifest_path, sep="\t")
                work = [(r["sample_id"], r["bam_path"]) for _, r in manifest.iterrows()]
                hap_rows = []
                with _cf.ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
                    futs = {ex.submit(haplotype_score, bp, pairs): sid for sid, bp in work}
                    for i, fut in enumerate(_cf.as_completed(futs), 1):
                        sid = futs[fut]
                        sc = fut.result()
                        sc["sample_id"] = sid
                        hap_rows.append(sc)
                        if i % 5 == 0 or i == len(work):
                            log(f"  haplotypes done {i}/{len(work)}")
                haplotypes_df = pd.DataFrame(hap_rows).set_index("sample_id")
            except Exception as exc:
                log(f"  Skipping haplotypes: {exc}")

            try:
                log("Computing per-sample molecule counts")
                scout = pd.read_csv(scout_sites_path, sep="\t")
                manifest = pd.read_csv(bam_manifest_path, sep="\t")
                work = [(r["sample_id"], r["bam_path"]) for _, r in manifest.iterrows()]
                mol_rows = []
                with _cf.ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
                    futs = {ex.submit(_mc_remote, sid, bp, scout, 500): sid for sid, bp in work}
                    for i, fut in enumerate(_cf.as_completed(futs), 1):
                        mc = fut.result()
                        mol_rows.append(mc)
                        if i % 5 == 0 or i == len(work):
                            log(f"  molecule_counts done {i}/{len(work)}")
                molecules_df = pd.DataFrame(mol_rows).set_index("sample_id")
            except Exception as exc:
                log(f"  Skipping molecule counts: {exc}")

    log("Assembling per-sample and per-pair score tables")
    per_sample = assemble_per_sample(
        samples=list(gt_m.columns),
        sample_to_group=sample_to_group,
        n_called_sites=n_called,
        directional=directional,
        charr=charr_df,
        tail=tail_df,
        haplotypes=haplotypes_df,
        molecules=molecules_df,
        umi_jaccard=umi_jaccard,
        umi_background=umi_background,
    )
    per_pair = assemble_per_pair(directional, sample_to_group, umi_jaccard)
    write_reports(per_sample, per_pair, out_dir=out_dir)

    directional["score"].to_csv(out_dir / "directional_contamination.tsv", sep="\t", float_format="%.6f")
    directional["score_het"].to_csv(out_dir / "directional_contamination_het_only.tsv", sep="\t", float_format="%.6f")
    directional["score_homalt"].to_csv(out_dir / "directional_contamination_homalt_only.tsv", sep="\t", float_format="%.6f")
    directional["n_informative_homalt"].to_csv(out_dir / "n_informative_sites_homalt.tsv", sep="\t")
    directional["background"].to_csv(out_dir / "background_alt_vaf.tsv", sep="\t", float_format="%.6f")
    ident.to_csv(out_dir / "identity_concordance_loh.tsv", sep="\t", float_format="%.4f")
    if not per_chrom_df.empty:
        per_chrom_df.to_csv(out_dir / "per_chrom_scores.tsv", sep="\t", float_format="%.6f", index=False)

    log(f"Wrote score tables to {out_dir}/")
    return {
        "per_sample": per_sample,
        "per_pair": per_pair,
        "directional": directional,
        "identity": ident,
        "sample_to_group": sample_to_group,
        "charr": charr_df,
        "tails": tail_df,
        "haplotypes": haplotypes_df,
        "molecules": molecules_df,
        "umi_jaccard": umi_jaccard,
        "per_chrom": per_chrom_df,
    }


def main():
    ap = argparse.ArgumentParser(description="Sentinel scoring stage")
    ap.add_argument("--ad-parquet", type=Path, default=WORK_DIR / "ad_long.parquet")
    ap.add_argument("--out-dir", type=Path, default=RESULTS_DIR)
    ap.add_argument("--min-depth-call", type=int, default=MIN_DEPTH_CALL)
    ap.add_argument("--min-depth-recipient", type=int, default=MIN_DEPTH_RECIPIENT)
    ap.add_argument("--het-lo", type=float, default=HET_CALL_LO)
    ap.add_argument("--het-hi", type=float, default=HET_CALL_HI)
    ap.add_argument("--hom-alt-lo", type=float, default=HOM_ALT_LO)
    ap.add_argument("--identity-threshold", type=float, default=IDENTITY_CLUSTER_THRESHOLD)
    ap.add_argument("--scout-sites", type=Path, default=CACHE_DIR / "scout_sites.tsv")
    ap.add_argument("--bam-manifest", type=Path, default=WORK_DIR / "run_manifest.tsv")
    ap.add_argument("--umi-length", type=int, default=10)
    ap.add_argument("--skip-bam-features", action="store_true")
    args = ap.parse_args()

    run(
        ad_parquet=args.ad_parquet,
        out_dir=args.out_dir,
        min_depth_call=args.min_depth_call,
        min_depth_recipient=args.min_depth_recipient,
        het_lo=args.het_lo, het_hi=args.het_hi, hom_alt_lo=args.hom_alt_lo,
        identity_threshold=args.identity_threshold,
        scout_sites_path=args.scout_sites,
        bam_manifest_path=args.bam_manifest,
        skip_bam_features=args.skip_bam_features,
        umi_length=args.umi_length,
    )


if __name__ == "__main__":
    main()

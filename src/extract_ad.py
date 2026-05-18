"""For each BAM, extract REF/ALT allele depths at the site list.

Per-sample parquet shards are written to <work>/ad_per_sample/. Downstream
loaders consume the shards directly as a Parquet dataset (no merge pass).

Parallelism: one worker per BAM (BAMs are large; sharing IO is what matters).
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from config import CACHE_DIR, MIN_BASE_QUAL, MIN_MAP_QUAL, WORK_DIR

BASES = "ACGT"
BASE_IDX = {b: i for i, b in enumerate(BASES)}

# Sites are clustered when consecutive positions are within MAX_GAP bp.
# Each cluster becomes one count_coverage fetch with margin MARGIN bp on either side.
MAX_GAP = 50_000
MARGIN = 200


def _cluster_positions(positions: np.ndarray):
    """Yield (start, stop, indices_into_positions) for each dense cluster of sites."""
    if positions.size == 0:
        return
    order = np.argsort(positions, kind="stable")
    sorted_pos = positions[order]
    gaps = np.diff(sorted_pos)
    break_points = np.where(gaps > MAX_GAP)[0] + 1
    cluster_slices = np.split(np.arange(sorted_pos.size), break_points)
    for sl in cluster_slices:
        if sl.size == 0:
            continue
        local_pos = sorted_pos[sl]
        idx_in_input = order[sl]
        yield int(local_pos.min()), int(local_pos.max()), local_pos, idx_in_input


def extract_one_bam(args):
    sample_id, bam_path, sites_path = args
    sites = pd.read_parquet(sites_path)
    out = np.zeros((len(sites), 4), dtype=np.int32)
    by_chrom = sites.groupby("chrom", sort=False)

    def read_ok(r):
        return (
            not r.is_unmapped and not r.is_secondary
            and not r.is_supplementary and not r.is_duplicate
            and not r.is_qcfail and r.mapping_quality >= MIN_MAP_QUAL
        )

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, sub in by_chrom:
            row_ix_all = sub.index.to_numpy()
            positions = sub["pos"].to_numpy()
            for cluster_min, cluster_max, local_pos, idx_in_chrom in _cluster_positions(positions):
                fetch_start = max(0, cluster_min - 1 - MARGIN)
                fetch_stop = cluster_max + MARGIN
                cov = bam.count_coverage(
                    contig=str(chrom),
                    start=fetch_start,
                    stop=fetch_stop,
                    quality_threshold=MIN_BASE_QUAL,
                    read_callback=read_ok,
                )
                counts = np.asarray(cov, dtype=np.int32)
                local_idx = local_pos - 1 - fetch_start
                picked = counts[:, local_idx].T
                out[row_ix_all[idx_in_chrom]] = picked

    rd = sites.copy()
    ref_idx = sites["ref"].map(BASE_IDX).to_numpy()
    alt_idx = sites["alt"].map(BASE_IDX).to_numpy()
    rd["ref_depth"] = out[np.arange(len(sites)), ref_idx]
    rd["alt_depth"] = out[np.arange(len(sites)), alt_idx]
    rd["other_depth"] = out.sum(axis=1) - rd["ref_depth"] - rd["alt_depth"]
    rd["depth"] = rd["ref_depth"] + rd["alt_depth"] + rd["other_depth"]
    rd.insert(0, "sample_id", sample_id)
    return rd[["sample_id", "chrom", "pos", "ref", "alt", "ref_depth", "alt_depth", "other_depth", "depth"]]


def prepare_sites_parquet(scout_tsv: Path, out_path: Path) -> Path:
    df = pd.read_csv(scout_tsv, sep="\t")
    keep = df[["chrom", "pos", "ref", "alt"]].copy()
    keep["chrom"] = keep["chrom"].astype(str)
    keep = keep.sort_values(["chrom", "pos"]).reset_index(drop=True)
    keep.to_parquet(out_path, index=False)
    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-proc", type=int, default=max(1, (os.cpu_count() or 4) // 2))
    parser.add_argument("--scout-tsv", type=Path, default=CACHE_DIR / "scout_sites.tsv")
    args = parser.parse_args()

    sites_path = CACHE_DIR / "sites.parquet"
    prepare_sites_parquet(args.scout_tsv, sites_path)
    n_sites = len(pd.read_parquet(sites_path))
    print(f"Loaded {n_sites:,} sites")

    manifest = pd.read_csv(WORK_DIR / "run_manifest.tsv", sep="\t")
    out_dir = WORK_DIR / "ad_per_sample"
    out_dir.mkdir(exist_ok=True)
    done = {p.stem for p in out_dir.glob("*.parquet")}
    work = [
        (r["sample_id"], r["bam_path"], str(sites_path))
        for _, r in manifest.iterrows()
        if r["sample_id"] not in done
    ]
    print(f"Already done: {len(done)} | To do: {len(work)} | Workers: {args.n_proc}")

    with mp.Pool(args.n_proc) as pool:
        for i, df in enumerate(pool.imap_unordered(extract_one_bam, work), 1):
            sid = df["sample_id"].iloc[0]
            df.to_parquet(out_dir / f"{sid}.parquet", index=False)
            mean_depth = df["depth"].mean()
            covered_pct = 100 * (df["depth"] >= 20).mean()
            print(f"[{i:>2}/{len(work)}] {sid}: mean_depth={mean_depth:.1f}  pct_d>=20={covered_pct:.1f}%")

    total_samples = len({p.stem for p in out_dir.glob("*.parquet")})
    print(f"\n{total_samples} per-sample parquet shards in {out_dir} ({n_sites:,} sites each)")


if __name__ == "__main__":
    main()

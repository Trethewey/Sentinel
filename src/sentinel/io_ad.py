"""Load the per-sample allele-depth parquet shards into a long-form DataFrame.

Reads either the per-sample shard directory written by `extract_ad.py`
(`<work>/ad_per_sample/`), or a single consolidated `<work>/ad_long.parquet`
file if the shard directory is absent.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
import pyarrow.dataset as ds

from .config import WORK_DIR


REQUIRED_COLS = ["sample_id", "chrom", "pos", "ref", "alt", "ref_depth", "alt_depth", "depth"]


def load_ad_long(
    path: Path = WORK_DIR / "ad_long.parquet",
    samples: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Read AD data. Prefers per-sample shard directory, falls back to a single parquet."""
    shard_dir = WORK_DIR / "ad_per_sample"
    if shard_dir.exists() and any(shard_dir.glob("*.parquet")):
        df = _read_shard_dir(shard_dir, samples)
    else:
        df = pd.read_parquet(path)
        if samples is not None:
            df = df[df["sample_id"].isin(list(samples))].reset_index(drop=True)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"AD table is missing columns: {missing}")
    df["chrom"] = df["chrom"].astype(str)
    return df


def load_ad_per_sample(
    directory: Path = WORK_DIR / "ad_per_sample",
    samples: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Read per-sample shards directly."""
    return _read_shard_dir(directory, samples)


def _read_shard_dir(directory: Path, samples: Optional[Iterable[str]]) -> pd.DataFrame:
    if samples is not None:
        wanted = set(samples)
        paths = [str(p) for p in sorted(directory.glob("*.parquet")) if p.stem in wanted]
    else:
        paths = [str(p) for p in sorted(directory.glob("*.parquet"))]
    if not paths:
        raise FileNotFoundError(f"No per-sample parquet shards found under {directory}")
    dataset = ds.dataset(paths, format="parquet")
    df = dataset.to_table().to_pandas()
    df["chrom"] = df["chrom"].astype(str)
    return df

"""UMI-derived metrics: per-sample molecule counts (F6) + cross-sample UMI footprint (F7).

F6 - molecule_counts
    Sentieon `--umi_tag XR --consensus` collapses raw reads into consensus
    reads. Consensus reads are the non-duplicate entries (their constituent
    raw reads keep FLAG 1024 = duplicate). Counting non-duplicate reads at
    scout positions therefore gives a molecule count, not a raw read count.
    Exposed per-sample to confirm the pipeline counts molecules and to give
    a downstream PCR duplication factor (raw / consensus).

F7 - umi_footprint_overlap
    For each pair of samples in a run, compute the Jaccard overlap of their
    UMI sequence sets (tag XR by default). Background overlap is ~ 4^-L for
    UMI length L; any pair above background by >= an order of magnitude
    suggests a sample swap or post-library contamination. Cheap to compute
    if we reservoir-sample reads per sample.
"""
from __future__ import annotations

import random
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import pysam
except ImportError:  # pragma: no cover - tests use synthetic UMI sets
    pysam = None

from config import UMI_MAX_PER_SAMPLE, UMI_TAG


# ---------------------------------------------------------------------------
# F6: per-sample molecule counts
# ---------------------------------------------------------------------------

def molecule_counts(
    bam_path: str,
    scout_sites: pd.DataFrame,
    sample_n_sites: Optional[int] = None,
    rng_seed: int = 0,
) -> dict:
    """Per-sample molecule-count summary at (optionally a sample of) scout sites.

    Counts consensus reads (non-duplicate, non-secondary, non-supplementary,
    non-qcfail) covering each scout site via pysam.AlignmentFile.count(),
    plus a raw-read tally (everything except secondary/supplementary/qcfail)
    so the caller can derive dup_factor = raw / consensus.
    """
    if pysam is None:
        raise RuntimeError("pysam not available")
    sites = scout_sites[["chrom", "pos"]].copy()
    sites["chrom"] = sites["chrom"].astype(str)
    if sample_n_sites and sample_n_sites < len(sites):
        sites = sites.sample(n=sample_n_sites, random_state=rng_seed).reset_index(drop=True)

    consensus = []
    raw = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for _, row in sites.iterrows():
            c = str(row["chrom"]); p = int(row["pos"])
            n_cons = bam.count(
                contig=c, start=p - 1, stop=p,
                read_callback=lambda r: (
                    not r.is_unmapped and not r.is_secondary and not r.is_supplementary
                    and not r.is_duplicate and not r.is_qcfail
                ),
            )
            n_raw = bam.count(
                contig=c, start=p - 1, stop=p,
                read_callback=lambda r: (
                    not r.is_unmapped and not r.is_secondary and not r.is_supplementary
                    and not r.is_qcfail
                ),
            )
            consensus.append(n_cons); raw.append(n_raw)
    cons_arr = np.array(consensus, dtype=np.int64)
    raw_arr = np.array(raw, dtype=np.int64)
    nonzero = cons_arr > 0
    dup_factor = float(raw_arr.sum() / cons_arr.sum()) if cons_arr.sum() > 0 else np.nan
    return {
        "n_sites_evaluated":        int(len(cons_arr)),
        "n_consensus_reads":        int(cons_arr.sum()),
        "n_raw_reads":              int(raw_arr.sum()),
        "mean_molecules_per_site":  float(cons_arr[nonzero].mean()) if nonzero.any() else 0.0,
        "median_molecules_per_site": float(np.median(cons_arr[nonzero])) if nonzero.any() else 0.0,
        "dup_factor":               dup_factor,
    }


# ---------------------------------------------------------------------------
# F7: cross-sample UMI footprint
# ---------------------------------------------------------------------------

def _read_umi(read, tag: str = UMI_TAG) -> Optional[str]:
    try:
        v = read.get_tag(tag)
    except KeyError:
        return None
    return str(v)


def collect_umis(
    bam_path: str,
    tag: str = UMI_TAG,
    max_umis: int = UMI_MAX_PER_SAMPLE,
    region: Optional[Tuple[str, int, int]] = None,
    rng_seed: int = 0,
) -> set:
    """Reservoir-sample distinct UMI strings from a BAM (set, not multiset)."""
    if pysam is None:
        raise RuntimeError("pysam not available")
    rng = random.Random(rng_seed)
    umis: List[str] = []
    seen = set()
    n_seen = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        it = bam.fetch(*region) if region else bam.fetch(until_eof=True)
        for read in it:
            if read.is_secondary or read.is_supplementary:
                continue
            u = _read_umi(read, tag)
            if u is None:
                continue
            if u in seen:
                continue
            seen.add(u)
            n_seen += 1
            if len(umis) < max_umis:
                umis.append(u)
            else:
                # reservoir replacement (Vitter's algorithm R)
                k = rng.randrange(n_seen)
                if k < max_umis:
                    umis[k] = u
    return set(umis)


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union > 0 else 0.0


def umi_footprint_overlap(
    umi_sets: Dict[str, set],
    sample_subset: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Pairwise UMI-set Jaccard. Accepts pre-collected umi_sets to keep the
    BAM I/O out of this function (testable + parallelisable).
    """
    samples = list(sample_subset) if sample_subset is not None else sorted(umi_sets.keys())
    n = len(samples)
    out = np.full((n, n), np.nan)
    for i in range(n):
        ai = umi_sets[samples[i]]
        out[i, i] = 1.0 if ai else np.nan
        for j in range(i + 1, n):
            v = jaccard(ai, umi_sets[samples[j]])
            out[i, j] = out[j, i] = v
    return pd.DataFrame(out, index=samples, columns=samples)


def background_jaccard_estimate(umi_length: int) -> float:
    """Expected baseline Jaccard between two independent UMI sets ~ 4^-L."""
    return float(4.0 ** (-umi_length))


def umi_footprint_from_bams(
    bam_paths_by_sample: Dict[str, str],
    tag: str = UMI_TAG,
    max_umis_per_sample: int = UMI_MAX_PER_SAMPLE,
    region: Optional[Tuple[str, int, int]] = None,
) -> Tuple[pd.DataFrame, Dict[str, set]]:
    """Convenience wrapper: collect UMI sets then compute Jaccard matrix.

    Returns (jaccard_matrix, umi_sets) so callers can reuse the sets.
    """
    sets: Dict[str, set] = {}
    for sid, bp in bam_paths_by_sample.items():
        sets[sid] = collect_umis(bp, tag=tag, max_umis=max_umis_per_sample, region=region)
    return umi_footprint_overlap(sets), sets

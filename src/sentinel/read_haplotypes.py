"""Read-level 3-haplotype counter on consensus reads (P1c + UMI angle #3).

For every pair of scout SNPs within HAP_PAIR_MAX_DIST bp of each other on the
same chromosome, examine consensus reads (and their mates) that span BOTH
positions and tally the observed (allele_a, allele_b) haplotypes.

A pure diploid sample emits at most 2 haplotypes per locus-pair. Three or
more haplotypes with non-trivial support is direct evidence of mixed source
DNA - and the signal does NOT require a donor to be on the run, so it
catches off-flowcell contamination too.

Consensus filtering rule (Sentieon `--consensus`): SKIP `is_duplicate`
because the consensus read is the surviving non-duplicate after UMI collapse.
Each consensus read = one independent original molecule, which makes the
3-haplotype count robust against PCR-amplified false support.

Public surface:
    find_site_pairs(scout_sites, max_dist=300) -> list[(site_a, site_b)]
    count_haplotypes_for_pair(bam_path, site_a, site_b, min_base_q=20) -> Counter
    haplotype_score(bam_path, site_pairs, max_pairs=2000) -> dict
"""
from __future__ import annotations

from collections import Counter
from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd

try:
    import pysam  # only required for BAM I/O
except ImportError:  # pragma: no cover - tests use a synthetic counter path
    pysam = None

from .config import HAP_MAX_PAIRS_EVAL, HAP_MIN_BASE_QUAL, HAP_MIN_SUPPORT_PER_HAP, HAP_PAIR_MAX_DIST

Site = Tuple[str, int, str, str]  # (chrom, pos_1based, ref, alt)


# ---------------------------------------------------------------------------
# Site-pair enumeration
# ---------------------------------------------------------------------------

def find_site_pairs(
    scout_sites: pd.DataFrame,
    max_dist: int = HAP_PAIR_MAX_DIST,
) -> List[Tuple[Site, Site]]:
    """Enumerate all ordered site pairs on the same chromosome within max_dist.

    Pairs are returned with pos_a < pos_b. scout_sites must have at minimum
    columns chrom, pos, ref, alt. Positions are 1-based.
    """
    df = scout_sites[["chrom", "pos", "ref", "alt"]].copy()
    df["chrom"] = df["chrom"].astype(str)
    df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)
    pairs: List[Tuple[Site, Site]] = []
    for _, sub in df.groupby("chrom", sort=False):
        sub = sub.reset_index(drop=True)
        positions = sub["pos"].to_numpy()
        n = len(sub)
        for i in range(n):
            pi = int(positions[i])
            # walk forward until distance exceeds max_dist
            for j in range(i + 1, n):
                pj = int(positions[j])
                if pj - pi > max_dist:
                    break
                if pj == pi:
                    continue
                a = (sub.at[i, "chrom"], pi, sub.at[i, "ref"], sub.at[i, "alt"])
                b = (sub.at[j, "chrom"], pj, sub.at[j, "ref"], sub.at[j, "alt"])
                pairs.append((a, b))
    return pairs


# ---------------------------------------------------------------------------
# Per-pair haplotype counting (BAM I/O)
# ---------------------------------------------------------------------------

def _base_at_position(read, ref_pos_0based: int, min_q: int):
    """Return the base call (A/C/G/T) at ref position, or None."""
    if read.query_sequence is None or read.query_qualities is None:
        return None
    # get_aligned_pairs(matches_only=True) yields (qpos, rpos) for aligned only
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if rpos == ref_pos_0based:
            if qpos is None:
                return None
            if read.query_qualities[qpos] < min_q:
                return None
            base = read.query_sequence[qpos]
            if base in "ACGT":
                return base
            return None
    return None


def count_haplotypes_for_pair(
    bam_path: str,
    site_a: Site,
    site_b: Site,
    min_base_q: int = HAP_MIN_BASE_QUAL,
):
    """Count (allele_a, allele_b) combinations across consensus reads / pairs.

    A fragment contributes if EITHER read in the pair (or a single spanning
    read) yields a base at BOTH positions. When the two positions are on
    different reads of the same pair, the bases are combined per fragment_id.
    """
    if pysam is None:
        raise RuntimeError("pysam not available")
    chrom_a, pos_a, _, _ = site_a
    chrom_b, pos_b, _, _ = site_b
    if chrom_a != chrom_b:
        raise ValueError("site pair must share chromosome")
    lo, hi = min(pos_a, pos_b), max(pos_a, pos_b)
    rpos_a = pos_a - 1
    rpos_b = pos_b - 1

    # qname -> {pos_a_base, pos_b_base}
    by_frag = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom_a, lo - 1, hi):
            if read.is_duplicate or read.is_secondary or read.is_supplementary:
                continue
            if read.is_unmapped or read.is_qcfail:
                continue
            qn = read.query_name
            frag = by_frag.setdefault(qn, {})
            # Only inspect a position if the read actually spans it
            if read.reference_start <= rpos_a < read.reference_end and "a" not in frag:
                b = _base_at_position(read, rpos_a, min_base_q)
                if b is not None:
                    frag["a"] = b
            if read.reference_start <= rpos_b < read.reference_end and "b" not in frag:
                b = _base_at_position(read, rpos_b, min_base_q)
                if b is not None:
                    frag["b"] = b

    counter: Counter = Counter()
    for frag in by_frag.values():
        if "a" in frag and "b" in frag:
            counter[(frag["a"], frag["b"])] += 1
    return counter


# ---------------------------------------------------------------------------
# Per-sample scoring (synthetic-counter-friendly)
# ---------------------------------------------------------------------------

def score_from_counters(
    counters: Iterable[Counter],
    min_support: int = HAP_MIN_SUPPORT_PER_HAP,
) -> dict:
    """Aggregate a sequence of per-pair haplotype Counters into a sample score.

    Drops haplotype combinations whose count < min_support (singleton noise),
    then for each remaining pair-with-evidence:
       - if >= 3 distinct surviving haplotypes -> count as flagged
       - minor_hap_fraction = sum(counts of non-top-2 haps) / total_supported_counts
    Returns dict with n_pairs_evaluated, n_pairs_with_3plus_haps, frac_3plus,
    mean_minor_hap_fraction (over flagged pairs only; NaN if none).
    """
    n_eval = 0
    n_flag = 0
    minor_fracs = []
    for c in counters:
        if not c:
            continue
        # keep only entries with support >= min_support
        kept = {h: n for h, n in c.items() if n >= min_support}
        if not kept:
            continue
        n_eval += 1
        n_haps = len(kept)
        if n_haps >= 3:
            n_flag += 1
            total = sum(kept.values())
            top2 = sum(sorted(kept.values(), reverse=True)[:2])
            minor = total - top2
            minor_fracs.append(minor / total if total > 0 else 0.0)
    frac_3plus = (n_flag / n_eval) if n_eval > 0 else np.nan
    mean_minor = float(np.mean(minor_fracs)) if minor_fracs else np.nan
    return {
        "n_pairs_evaluated": int(n_eval),
        "n_pairs_with_3plus_haps": int(n_flag),
        "frac_3plus": float(frac_3plus) if n_eval > 0 else np.nan,
        "mean_minor_hap_fraction": mean_minor,
    }


def haplotype_score(
    bam_path: str,
    site_pairs: List[Tuple[Site, Site]],
    max_pairs: int = HAP_MAX_PAIRS_EVAL,
    min_base_q: int = HAP_MIN_BASE_QUAL,
    min_support: int = HAP_MIN_SUPPORT_PER_HAP,
) -> dict:
    """End-to-end per-sample haplotype score from a BAM and a list of pairs."""
    pairs = site_pairs[:max_pairs] if (max_pairs and len(site_pairs) > max_pairs) else site_pairs
    counters = (
        count_haplotypes_for_pair(bam_path, a, b, min_base_q=min_base_q)
        for a, b in pairs
    )
    return score_from_counters(counters, min_support=min_support)

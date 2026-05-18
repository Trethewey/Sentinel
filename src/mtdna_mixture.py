"""mtDNA mixture detector.

Off-target chrM coverage on most targeted-capture panels is high enough
(1000-10000x) to detect mtDNA haplogroup mixing. Mitochondria are inherited maternally and
travel with the cell, so:

- Single individual (incl. SCT - all blood cells are donor's, single mtDNA)
  -> monoclonal mtDNA, almost no positions with mixed VAF
- Contamination (two cell populations mixed)
  -> mtDNA also mixed -> positions with VAF in (0.05, 0.95) above background

This is a single-sample assay - no cohort required. Output is the fraction of
high-depth chrM positions with mixed VAF, plus a binary "mixed?" call.

Used by sct_discriminator.py to disambiguate SCT vs contamination when the
nuclear directional matrix flags a sample.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pysam

# Possible chrM contig names across reference builds
CHR_M_CANDIDATES = ("chrM", "M", "MT", "chrMT", "Mt")
MIN_MT_DEPTH = 100
MIXED_VAF_LO = 0.05
MIXED_VAF_HI = 0.95
ERROR_RATE_BG = 0.005   # baseline single-base error rate per position
MIN_ALT_READS_FOR_MIXED = 5


def find_chrm_contig(bam: pysam.AlignmentFile) -> Optional[str]:
    refs = set(bam.references)
    for c in CHR_M_CANDIDATES:
        if c in refs:
            return c
    return None


def chrm_mixture_score(bam_path: Path, min_base_quality: int = 20,
                       min_mapping_quality: int = 20):
    """Walk chrM pileup; return dict with mixture stats."""
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        chrm = find_chrm_contig(bam)
        if chrm is None:
            return {"chrm": None, "n_positions": 0, "n_high_depth": 0,
                    "n_mixed": 0, "frac_mixed": np.nan, "mt_mean_depth": 0,
                    "verdict": "no_chrM"}
        cov = bam.count_coverage(
            contig=chrm,
            quality_threshold=min_base_quality,
            read_callback=lambda r: (
                not r.is_unmapped and not r.is_secondary
                and not r.is_supplementary and not r.is_duplicate
                and not r.is_qcfail and r.mapping_quality >= min_mapping_quality
            ),
        )
        counts = np.asarray(cov, dtype=np.int32)  # shape (4, length)
        depth = counts.sum(axis=0)
        n_positions = int((depth > 0).sum())
        if n_positions == 0:
            return {"chrm": chrm, "n_positions": 0, "n_high_depth": 0,
                    "n_mixed": 0, "frac_mixed": np.nan, "mt_mean_depth": 0,
                    "verdict": "no_mt_coverage"}
        ref = pysam.FastaFile  # not strictly needed; use the major allele as 'ref'
        major = counts.max(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            minor = depth - major
            minor_vaf = np.where(depth > 0, minor / np.where(depth > 0, depth, 1), 0)
        hd = depth >= MIN_MT_DEPTH
        n_high_depth = int(hd.sum())
        mixed_mask = hd & (minor_vaf >= MIXED_VAF_LO) & (minor_vaf <= MIXED_VAF_HI) & (minor >= MIN_ALT_READS_FOR_MIXED)
        n_mixed = int(mixed_mask.sum())
        frac_mixed = (n_mixed / n_high_depth) if n_high_depth else float("nan")
        mt_mean_depth = float(depth[depth > 0].mean())
        # Verdict
        if n_high_depth < 50:
            verdict = "insufficient_chrM_depth"
        elif frac_mixed > 0.005:
            verdict = "mtdna_mixed"
        else:
            verdict = "mtdna_monoclonal"
    return {
        "chrm": chrm,
        "n_positions": n_positions,
        "n_high_depth": n_high_depth,
        "n_mixed": n_mixed,
        "frac_mixed": frac_mixed,
        "mt_mean_depth": mt_mean_depth,
        "verdict": verdict,
    }

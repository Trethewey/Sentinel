"""Stem-cell-transplant vs contamination discriminator.

Combines three orthogonal signals to classify samples that look genetically
non-self into one of {SCT_LIKE, CONTAMINATION, INDETERMINATE, PASS}:

  1. Nuclear directional score (score_homalt)
     - High = strong non-self signal
     - Both SCT and heavy contamination score high here.
  2. Fraction of phaseable site pairs supporting 3+ haplotypes
     - SCT: near zero (donor is one diploid individual, at most 2 haplotypes
       per locus)
     - Contamination: positive (two diploid individuals mixed, 3-4 haplotypes
       at some loci)
  3. mtDNA mixture fraction
     - SCT: near zero (donor cells, single mtDNA haplogroup)
     - Contamination: positive (mixed cell populations, mixed mtDNA)

When the nuclear score is below threshold the sample is PASS. When it is
above threshold the other two signals discriminate SCT from contamination,
and disagreement between them marks the case as INDETERMINATE for review.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

# Default thresholds; tune on a labelled cohort if available.
NUCLEAR_SCORE_HIGH = 0.20
FRAC_3PLUS_HAP_LOW = 0.01
FRAC_3PLUS_HAP_HIGH = 0.05
MT_FRAC_MIXED_LOW = 0.005
MT_FRAC_MIXED_HIGH = 0.02

Verdict = Literal["PASS", "SCT_LIKE", "CONTAMINATION", "INDETERMINATE"]


@dataclass
class SCTCall:
    verdict: Verdict
    nuclear_score_homalt: float
    frac_3plus_hap: float
    mt_frac_mixed: float
    rationale: str


def call(nuclear_score_homalt: float,
         frac_3plus_hap: float,
         mt_frac_mixed: float) -> SCTCall:
    if nuclear_score_homalt < NUCLEAR_SCORE_HIGH:
        return SCTCall(
            verdict="PASS",
            nuclear_score_homalt=nuclear_score_homalt,
            frac_3plus_hap=frac_3plus_hap,
            mt_frac_mixed=mt_frac_mixed,
            rationale="nuclear non-self signal is below threshold; no further call",
        )

    hap_low = frac_3plus_hap < FRAC_3PLUS_HAP_LOW
    hap_high = frac_3plus_hap > FRAC_3PLUS_HAP_HIGH
    mt_low = mt_frac_mixed < MT_FRAC_MIXED_LOW
    mt_high = mt_frac_mixed > MT_FRAC_MIXED_HIGH

    if hap_low and mt_low:
        return SCTCall(
            verdict="SCT_LIKE",
            nuclear_score_homalt=nuclear_score_homalt,
            frac_3plus_hap=frac_3plus_hap,
            mt_frac_mixed=mt_frac_mixed,
            rationale="strong nuclear non-self + monoclonal haplotypes + monoclonal mtDNA "
                      "consistent with stem cell transplant (donor cells, not a mixture)",
        )
    if hap_high and mt_high:
        return SCTCall(
            verdict="CONTAMINATION",
            nuclear_score_homalt=nuclear_score_homalt,
            frac_3plus_hap=frac_3plus_hap,
            mt_frac_mixed=mt_frac_mixed,
            rationale="strong nuclear non-self + >2-haplotype evidence + mixed mtDNA "
                      "consistent with cross-sample contamination",
        )
    return SCTCall(
        verdict="INDETERMINATE",
        nuclear_score_homalt=nuclear_score_homalt,
        frac_3plus_hap=frac_3plus_hap,
        mt_frac_mixed=mt_frac_mixed,
        rationale=(f"signals disagree (3plus_hap={frac_3plus_hap:.3f}, "
                   f"mt_mixed={mt_frac_mixed:.4f}); manual review needed. "
                   "Possible explanations: tumour subclones, mtDNA heteroplasmy, "
                   "partial chimerism, or post-PCR contamination."),
    )

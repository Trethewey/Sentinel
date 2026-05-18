"""Project paths, constants, and tunable thresholds.

Paths default to the current working directory. The CLI (`sentinel run`)
overrides them via environment variables before calling each stage, so users
do not normally edit this file.
"""
from __future__ import annotations

import os
from pathlib import Path

# --- Paths -------------------------------------------------------------------
PROJECT_DIR = Path(os.environ.get("SENTINEL_PROJECT", "."))
BAMS_DIR    = Path(os.environ.get("SENTINEL_BAMS",   PROJECT_DIR / "bams"))
REF_FA      = Path(os.environ.get("SENTINEL_REF",    PROJECT_DIR / "reference.fa"))
WORK_DIR    = Path(os.environ.get("SENTINEL_WORK",   PROJECT_DIR / "work"))
CACHE_DIR   = Path(os.environ.get("SENTINEL_CACHE",  PROJECT_DIR / "cache"))
RESULTS_DIR = Path(os.environ.get("SENTINEL_RESULTS", PROJECT_DIR / "results"))

# Note: directories are created on demand by the stages that need them,
# not at import time, so importing config has no side effects on the filesystem.

# --- Chromosomes -------------------------------------------------------------
AUTOSOMES = [str(i) for i in range(1, 23)]
SEX_CHROMS = ["X", "Y"]
ALL_NUCLEAR = AUTOSOMES + SEX_CHROMS

# --- Read / base QC ----------------------------------------------------------
MIN_BASE_QUAL = 20
MIN_MAP_QUAL = 20

# --- Scout / genotyping ------------------------------------------------------
MIN_DEPTH_SCOUT = 40
MIN_DEPTH_CALL = 20
MIN_DEPTH_RECIPIENT = 20
HET_VAF_RANGE = (0.30, 0.70)
HET_CALL_LO = 0.25
HET_CALL_HI = 0.75
HOM_ALT_LO = 0.90
MAX_VAF_HOM_REF = 0.05

# --- Identity clustering -----------------------------------------------------
IDENTITY_CLUSTER_THRESHOLD = 0.95
LOH_ALT_EVIDENCE_MIN_READS = 1

# --- Tail-anomaly bands ------------------------------------------------------
TAIL_LO = (0.05, 0.15)
TAIL_HI = (0.85, 0.95)

# --- Read-haplotype pairs ----------------------------------------------------
HAP_PAIR_MAX_DIST = 300
HAP_MIN_BASE_QUAL = 20
HAP_MIN_SUPPORT_PER_HAP = 2
HAP_MAX_PAIRS_EVAL = 2000

# --- UMI metrics -------------------------------------------------------------
UMI_TAG = "XR"
UMI_MAX_PER_SAMPLE = 200_000

# --- Verdict thresholds (initial values; tunable) ---------------------------
VERDICT = {
    "fail_score_homalt":       0.01,
    "fail_n_informative_min":  50,
    "fail_homalt_deflation":   0.02,
    "fail_frac_3plus_haps":    0.05,
    "warn_score_homalt":       0.005,
    "warn_homalt_deflation":   0.01,
    "warn_frac_3plus_haps":    0.025,
    "warn_vaf_tail_fraction":  0.05,
    "warn_umi_overlap_mult":   10.0,
}

# --- Genotype constants ------------------------------------------------------
GT_NOCALL, GT_HOM_REF, GT_HET, GT_HOM_ALT = 0, 1, 2, 3
GT_LABEL = {0: "no_call", 1: "hom_ref", 2: "het", 3: "hom_alt"}

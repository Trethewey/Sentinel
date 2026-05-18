"""Sentinel CLI.

Subcommands:
  build-panel   build a SNP catalog by intersecting a panel BED with a master VCF
  run           extract allele depths, score, post-process, render the report
  report        re-render the HTML report from an existing results directory
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

from . import __version__


def _set_env(env: dict) -> dict:
    """Merge env over os.environ and apply to current process."""
    for k, v in env.items():
        os.environ[k] = str(v)
    return os.environ


def cmd_build_panel(args: argparse.Namespace) -> int:
    from . import panel_builder
    parser_args = argparse.Namespace(
        bed=str(args.bed),
        master_db=str(args.master_db),
        panel_name=args.panel_name,
        out=str(args.out),
    )
    panel_builder.from_bed(parser_args)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    work = out / "work"
    cache = out / "cache"
    results = out / "results"
    scores = results / "scores"
    for p in (work, cache, results, scores, work / "ad_per_sample"):
        p.mkdir(parents=True, exist_ok=True)

    catalog_dst = cache / "scout_sites.tsv"
    _stage_catalog(args.catalog, catalog_dst, args.strip_chr_prefix)

    manifest_tsv = work / "run_manifest.tsv"
    _build_manifest(args.sample_sheet, args.bams_dir, manifest_tsv)

    _set_env({
        "PROJ": str(out),
        "SENTINEL_PROJECT": str(out),
        "SENTINEL_WORK": str(work),
        "SENTINEL_CACHE": str(cache),
        "SENTINEL_RESULTS": str(results),
        "SENTINEL_REF": str(args.ref) if args.ref else "",
    })

    # Need to re-import config after env vars are set so paths pick them up.
    from . import config

    # Stage 1: extract allele depths from each BAM at the catalog sites.
    from .extract_ad import main as extract_main
    sys.argv = ["extract_ad", "--n-proc", str(args.threads), "--scout-tsv", str(catalog_dst)]
    extract_main()

    # Stage 2: score (identity, directional, deflation, tails, haplotypes).
    from .pipeline import run as pipeline_run
    pipeline_run(
        ad_parquet=work / "ad_long.parquet",
        out_dir=scores,
        scout_sites_path=catalog_dst,
        bam_manifest_path=manifest_tsv,
    )

    # Stage 3: post-process verdicts.
    from .post_process import main as post_main
    post_main()

    # Stage 4: HTML + XLSX reports.
    from .make_html_report import main as html_main
    html_main()
    from .make_xlsx_report import main as xlsx_main
    xlsx_main()
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    _set_env({"PROJ": str(Path(args.results_dir).resolve())})
    from .make_html_report import main as html_main
    from .make_xlsx_report import main as xlsx_main
    html_main()
    xlsx_main()
    return 0


def _stage_catalog(src: Path, dst: Path, strip_chr: bool) -> None:
    """Copy catalog TSV to cache, optionally stripping the leading 'chr' from contigs."""
    with open(src, encoding="utf-8") as f_in, open(dst, "w", encoding="utf-8", newline="") as f_out:
        next(f_in)  # header
        f_out.write("chrom\tpos\tref\talt\tn_scout_donors\tscout_donors\n")
        for line in f_in:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            if strip_chr and chrom.startswith("chr"):
                chrom = chrom[3:]
            f_out.write(f"{chrom}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t1\tpanel\n")


def _build_manifest(sample_sheet: Path, bams_dir: Path, out_tsv: Path) -> None:
    rows = []
    with open(sample_sheet, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            sid = r["sample_id"].strip()
            if not sid:
                continue
            bam = _find_bam(bams_dir, sid)
            rows.append({
                "sample_id": sid,
                "bam_path": str(bam),
                "lab_id": sid,
                "patient_id": (r.get("patient_id") or "").strip() or sid,
                "sample_type": (r.get("sample_type") or "").strip() or "clinical",
                "well": (r.get("well") or "").strip(),
                "expected_identity": (r.get("expected_identity") or "").strip() or sid,
                "replicate": (r.get("replicate") or "").strip(),
                "sex_declared": (r.get("sex_declared") or "U").strip(),
            })
    cols = ["sample_id", "bam_path", "lab_id", "patient_id", "sample_type",
            "well", "expected_identity", "replicate", "sex_declared"]
    with open(out_tsv, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(r[c] for c in cols) + "\n")


def _find_bam(bams_dir: Path, sample_id: str) -> Path:
    for suffix in ("_sorted.bam", ".bam"):
        p = bams_dir / f"{sample_id}{suffix}"
        if p.exists():
            return p
    raise FileNotFoundError(f"No BAM found for {sample_id} under {bams_dir}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="sentinel",
                                 description="Cross-sample contamination and identity detector for targeted NGS runs")
    ap.add_argument("--version", action="version", version=f"sentinel {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    pb = sub.add_parser("build-panel", help="Build a SNP catalog from a panel BED")
    pb.add_argument("--bed", type=Path, required=True, help="Panel target BED file")
    pb.add_argument("--master-db", type=Path, required=True,
                    help="Master SNP VCF (bgzipped + tabix-indexed, e.g. gnomAD common SNPs)")
    pb.add_argument("--panel-name", required=True, help="Panel name (used in output filenames)")
    pb.add_argument("--out", type=Path, required=True, help="Output directory for catalog + metadata")
    pb.set_defaults(func=cmd_build_panel)

    rn = sub.add_parser("run", help="Run the detector against a sample sheet")
    rn.add_argument("--sample-sheet", type=Path, required=True,
                    help="CSV with at least a sample_id column")
    rn.add_argument("--bams-dir", type=Path, required=True,
                    help="Directory containing <sample_id>_sorted.bam files")
    rn.add_argument("--catalog", type=Path, required=True,
                    help="Site catalog TSV from build-panel")
    rn.add_argument("--ref", type=Path, default=None,
                    help="Reference FASTA (optional)")
    rn.add_argument("--out", type=Path, required=True, help="Output directory")
    rn.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 4) // 2))
    rn.add_argument("--strip-chr-prefix", action="store_true",
                    help="Strip the leading 'chr' from catalog contigs to match no-prefix BAMs")
    rn.set_defaults(func=cmd_run)

    rp = sub.add_parser("report", help="Re-render the HTML report from an existing results directory")
    rp.add_argument("--results-dir", type=Path, required=True)
    rp.set_defaults(func=cmd_report)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

"""Panel builder: build a Sentinel SNP catalog by intersecting a panel BED
with a master SNP database VCF.

Outputs per panel:
  <out_dir>/<panel_name>_sites.tsv      tab-separated chrom / pos / ref / alt / af
  <out_dir>/<panel_name>_meta.json      panel name, build params, software version, counts

Usage:
  panel_builder.py from-bed \\
      --bed panel.bed \\
      --master-db /path/to/master_snps.vcf.gz \\
      --panel-name MyPanel \\
      --out /path/to/catalogs/
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import subprocess
from pathlib import Path


SENTINEL_VERSION = "0.1.0"


def md5_of(path: str, chunk: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def detect_db_contig_style(db_path: str) -> str:
    """Inspect the VCF header to see if contigs are 'chr1' or '1'."""
    cmd = ["bcftools", "view", "-h", db_path]
    hdr = subprocess.check_output(cmd, text=True)
    for line in hdr.splitlines():
        if line.startswith("##contig="):
            if "ID=chr" in line:
                return "chr"
            return "nochr"
    return "unknown"


def detect_bed_contig_style(bed_path: str) -> str:
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            first = line.split("\t")[0]
            return "chr" if first.startswith("chr") else "nochr"
    return "unknown"


def normalise_bed(bed_path: str, out_path: str, target_style: str) -> None:
    """Rewrite a BED with contigs matching target_style ('chr' or 'nochr')."""
    with open(bed_path) as fi, open(out_path, "w") as fo:
        for line in fi:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            ch = parts[0]
            if target_style == "chr" and not ch.startswith("chr"):
                ch = "chr" + ch
            elif target_style == "nochr" and ch.startswith("chr"):
                ch = ch[3:]
            parts[0] = ch
            fo.write("\t".join(parts) + "\n")


def from_bed(args):
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    sites_tsv = out_dir / f"{args.panel_name}_sites.tsv"
    meta_json = out_dir / f"{args.panel_name}_meta.json"

    print(f"[panel_builder] panel={args.panel_name}", flush=True)
    print(f"[panel_builder] master DB: {args.master_db}", flush=True)
    print(f"[panel_builder] panel BED: {args.bed}", flush=True)

    db_style = detect_db_contig_style(args.master_db)
    bed_style = detect_bed_contig_style(args.bed)
    print(f"[panel_builder] DB uses '{db_style}' contigs, BED uses '{bed_style}'", flush=True)

    tmp_bed = out_dir / f".{args.panel_name}_bed.bed"
    if db_style != bed_style:
        print(f"[panel_builder] normalising BED contigs to '{db_style}' style", flush=True)
        normalise_bed(args.bed, str(tmp_bed), db_style)
        use_bed = str(tmp_bed)
    else:
        use_bed = args.bed

    print("[panel_builder] intersecting master DB with panel BED", flush=True)
    cmd_view = ["bcftools", "view", "-R", use_bed, "-O", "v", args.master_db]
    cmd_query = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n"]
    n_sites = 0
    with open(sites_tsv, "w") as fo:
        fo.write("chrom\tpos\tref\talt\taf\n")
        p1 = subprocess.Popen(cmd_view, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd_query, stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
        p1.stdout.close()
        for line in p2.stdout:
            fo.write(line)
            n_sites += 1
        p2.wait()
        p1.wait()
    print(f"[panel_builder] wrote {sites_tsv}: {n_sites:,} sites", flush=True)

    if tmp_bed.exists():
        tmp_bed.unlink()

    per_chrom = {}
    with open(sites_tsv) as f:
        next(f)
        for line in f:
            ch = line.split("\t", 1)[0]
            per_chrom[ch] = per_chrom.get(ch, 0) + 1

    meta = {
        "panel_name": args.panel_name,
        "software": "Sentinel",
        "software_version": SENTINEL_VERSION,
        "build_date": datetime.datetime.utcnow().isoformat() + "Z",
        "build_mode": "from-bed",
        "master_db": args.master_db,
        "master_db_md5": md5_of(args.master_db) if Path(args.master_db).stat().st_size < 5 * 1024**3 else "skipped_large_file",
        "panel_bed": args.bed,
        "panel_bed_md5": md5_of(args.bed),
        "n_sites_total": n_sites,
        "n_sites_per_chrom": per_chrom,
        "min_maf": "inherited from master DB",
    }
    with open(meta_json, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[panel_builder] wrote {meta_json}", flush=True)
    print(f"[panel_builder] DONE  panel={args.panel_name}  sites={n_sites:,}", flush=True)


def main():
    parser = argparse.ArgumentParser(description="Build a Sentinel panel catalog.")
    subs = parser.add_subparsers(dest="cmd", required=True)

    p_bed = subs.add_parser("from-bed", help="Intersect a panel BED with the master SNP database")
    p_bed.add_argument("--bed", required=True, help="Panel target BED file")
    p_bed.add_argument("--master-db", required=True, help="Master SNP database VCF (bgzipped + indexed)")
    p_bed.add_argument("--panel-name", required=True, help="Panel name (used in output filenames)")
    p_bed.add_argument("--out", required=True, help="Output directory")
    p_bed.set_defaults(func=from_bed)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

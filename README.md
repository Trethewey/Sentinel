<p align="center">
  <img src="assets/logo.svg" alt="Sentinel" width="540">
</p>

Sentinel finds cross-sample contamination, sample swaps, and identity drift on
a run of targeted NGS samples. It works by comparing each sample's genotype
calls to every other sample on the same run, looking for the asymmetric
read-mixing signature that contamination leaves behind.

It is a single command that takes a BAM directory, a sample sheet, and a SNP
catalog, and writes a self-contained HTML report. No database, no daemon, no
sample-sheet-format gymnastics.

---

## Install

The detector ships as a Python package with a small CLI. Pick one:

```bash
# pip
pip install sentinel-cd

# conda (bioconda)
conda install -c bioconda sentinel-cd

# Docker
docker pull ghcr.io/<org>/sentinel:latest
```

`samtools` and `bcftools` must be on `PATH` (the conda and Docker images bundle
them).

Tested with Python 3.10-3.12 on Linux and WSL. macOS works but isn't part of
the CI matrix.

---

## What you need

| | |
|---|---|
| Aligned BAMs | one per sample, indexed (`.bai`), GRCh38 reference |
| Sample sheet | CSV, one row per sample, at minimum a `sample_id` column |
| Site catalog | TSV from `sentinel build-panel`, scoped to your panel |
| Reference FASTA | GRCh38, only required if you run the haplotype branch |

The site catalog is a one-off per panel: build it once, reuse it for every
run that uses the same baits.

---

## Quick start

```bash
# 1. Build the catalog for your panel (one-off)
sentinel build-panel \
    --bed mypanel.bed \
    --master-db gnomad_v41_common.vcf.gz \
    --panel-name MyPanel \
    --out catalogs/

# 2. Run the detector on a run
sentinel run \
    --sample-sheet samples.csv \
    --bams-dir /path/to/bams \
    --catalog catalogs/MyPanel_sites.tsv \
    --ref GRCh38.fa \
    --out run_2026_05_17/
```

After `sentinel run` finishes you'll find in `run_2026_05_17/results/`:

- `report.html` - the interactive single-file report (open in a browser)
- `report.xlsx` - the same per-sample table, plus a column legend, as Excel
- `per_sample_report.tsv` - one row per sample, machine-readable
- `per_pair_report.tsv` - long-form donor-to-recipient table for flagged pairs

Open `report.html` and you get a 96-well plate view, a flagged-samples table,
a per-sample contamination bar plot, and a sortable all-samples table. Click
any row to isolate that sample on the plate.

---

## Sample sheet

Minimum: one column.

```csv
sample_id
SAMPLE001
SAMPLE002
SAMPLE003
```

BAMs are matched in `--bams-dir` by `<sample_id>_sorted.bam` or
`<sample_id>.bam`.

Optional columns:

| Column | Purpose |
|---|---|
| `sample_type` | One of `control`, `reference_fresh`, `reference_ffpe`, `mix`, `clinical`. Drives the verdict thresholds; controls are held to a different bar than clinical samples. Defaults to `clinical` if absent. |
| `well` | 96-well coordinate (e.g. `A1`, `H12`). When present, the report shows a plate-layout view that highlights neighbouring-well cross-contamination. |
| `expected_identity` | The sample identity each BAM is labelled as. If the genetics don't match, the sample is flagged as a swap. Defaults to `sample_id`. |
| `patient_id` | Trio / family grouping. Only used if multiple samples share a patient_id. |

A template is in `examples/sample_sheet_template.csv` and an extended example
with the optional columns is in `examples/sample_sheet_with_optional.csv`.

---

## Commands

| Command | Purpose |
|---|---|
| `sentinel build-panel` | Build a SNP catalog by intersecting your panel BED with a master VCF |
| `sentinel run` | Extract allele depths, score, post-process, render reports |
| `sentinel report` | Re-render reports from an existing results directory |

Each command has `--help` with the full set of options.

---

## The master SNP database

`sentinel build-panel` needs a master VCF of common SNPs. The recommended
source is gnomAD v4.1, filtered to biallelic SNPs at MAF >= 5%.

Two pre-built distributions are published:

| File | Source | Sites | Size | Host |
|---|---|---|---|---|
| `gnomad_v41_exomes_common.vcf.gz` | gnomAD v4.1 exomes, MAF >= 5%, biallelic SNPs | ~2.6M | ~1.8 GB | GitHub Releases |
| `gnomad_v41_genomes_maf01.vcf.gz` | gnomAD v4.1 genomes, MAF >= 1%, biallelic SNPs | ~15M | ~15 GB | Zenodo, DOI [10.5281/zenodo.20264481](https://doi.org/10.5281/zenodo.20264481) |

The exomes file works for most capture panels (gene panels, exome panels,
targeted hybridisation-capture assays). Use the genomes file when your panel
includes non-coding bait regions (custom panels with intronic content, WGS,
etc.).

Both files are bgzipped and tabix-indexed. Download, drop next to your
analysis, point `sentinel build-panel --master-db` at them.

```bash
# exomes (1.8 GB) + index, from GitHub
curl -LO https://github.com/Trethewey/Sentinel/releases/download/v0.1.0/gnomad_v41_exomes_common.vcf.gz
curl -LO https://github.com/Trethewey/Sentinel/releases/download/v0.1.0/gnomad_v41_exomes_common.vcf.gz.tbi

# genomes (15 GB) + index, from Zenodo
curl -LO https://zenodo.org/records/20264481/files/gnomad_v41_genomes_maf01.vcf.gz
curl -LO https://zenodo.org/records/20264481/files/gnomad_v41_genomes_maf01.vcf.gz.tbi
```

To build your own from scratch see `tools/build_master_db.sh`.

---

## Docker

```bash
docker run --rm -v $PWD:/data ghcr.io/<org>/sentinel:latest run \
    --sample-sheet /data/samples.csv \
    --bams-dir /data/bams \
    --catalog /data/catalogs/MyPanel_sites.tsv \
    --ref /data/GRCh38.fa \
    --out /data/run_output
```

The image is based on `python:3.11-slim` with `samtools`, `bcftools`, and all
Python dependencies preinstalled.

---

## DNA Nexus

A DNA Nexus app wrapper is shipped in `dxapp/`. Build with:

```bash
cd dxapp && dx build --upload
```

The app exposes the same inputs as the CLI: sample sheet, BAMs, catalog,
reference FASTA, and the optional thread count. Outputs land back in the
DNA Nexus project as `report.html`, `per_sample_report.csv`, and
`per_pair_report.csv`.

---

## Output schema

`per_sample_report.tsv` columns (the headline ones):

| Column | Meaning |
|---|---|
| `sample_id` | From the sample sheet |
| `well` | Plate well (if supplied) |
| `sample_type` | From the sample sheet (defaults to `clinical`) |
| `verdict` | `PASS`, `WARN`, or `FAIL` |
| `top_score_homalt` | Headline contamination score; approximately equals the contamination fraction |
| `top_donor_sample_id` | Best guess at the source of the contamination |
| `homalt_deflation` | Drop in alt-VAF at the sample's own hom-alt sites |
| `identity_match` | `TRUE` if the genetics match the declared identity, `FALSE` if swapped |
| `concordance_to_expected` | Genetic similarity to the declared identity (0-1) |

The Excel report carries every column and a legend sheet explaining each one.

---

## How it works

At each SNP site in the catalog, Sentinel asks: at sites where the candidate
donor has the alt allele but the recipient should be hom-ref, how often does
the recipient actually carry alt reads? Pure samples score near zero;
contaminated samples carry the donor's alt reads at a fraction approximately
equal to the contamination level (alpha).

The detector layers this asymmetric matrix score with a CHARR-style hom-alt
deflation metric, a VAF-tail anomaly score, a 3+-haplotype read-level check,
and an identity-anchor cross-check. The verdict logic combines all of these
with sample-type-aware thresholds: controls and FFPE samples are held to
looser bars than fresh clinical samples.

For more detail, see `docs/method.md`.

---

## Citation

If you use Sentinel in published work, please cite the GitHub repository
(a preprint is in preparation) and the papers that introduced the methods
Sentinel builds on.

**Methods**

- Lu W, Gauthier LD, Poterba T et al. *CHARR efficiently estimates
  contamination from DNA sequencing samples.* Bioinformatics, 2023.
  [`charr.py`]
- Jun G, Flickinger M, Hetrick KN et al. *Detecting and estimating
  contamination of human DNA samples in sequencing and array-based genotype
  data.* Am J Hum Genet 91:839-848, 2012. [conceptual baseline for
  per-sample contamination QC]
- Pedersen BS, Bhetariya PJ, Brown J et al. *Somalier: rapid relatedness
  estimation for cancer and germline studies using efficient genome
  sketches.* Genome Medicine 12:62, 2020. [`identity_loh.py`,
  `anchor_check.py`]
- Schmitt MW, Kennedy SR, Salk JJ et al. *Detection of ultra-rare mutations
  by next-generation sequencing.* PNAS 109:14508-14513, 2012.
  [`read_haplotypes.py`]
- Weissensteiner H, Forer L, Fendt L et al. *Contamination detection in
  sequencing studies using the mitochondrial phylogeny (Haplocheck).*
  Genome Research 31:309-316, 2021. [`mtdna_mixture.py`,
  `sct_discriminator.py`]

**Reference data**

- Chen S, Francioli LC, Goodrich JK et al. *A genomic mutational constraint
  map using variation in 76,156 human genomes* (gnomAD v4). Nature
  625:92-100, 2024.
- Trethewey C. *Sentinel reference panel: gnomAD v4.1 genomes biallelic SNPs
  at MAF >= 1%* (2026). Zenodo. DOI: 10.5281/zenodo.20264481.

**Tools**

- Li H. *The Sequence Alignment/Map format and SAMtools.* Bioinformatics
  25:2078-2079, 2009.
- Danecek P, Bonfield JK, Liddle J et al. *Twelve years of SAMtools and
  BCFtools.* GigaScience 10:giab008, 2021.
- pysam: github.com/pysam-developers/pysam

**Scientific Python stack** (if your venue requires it): NumPy (Harris et
al. *Nature* 585:357-362, 2020), pandas (McKinney, *SciPy* 2010),
PyArrow, Plotly.

---

## License

See `LICENSE`.

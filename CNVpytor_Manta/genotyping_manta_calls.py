#!/usr/bin/env python3
# Script to genotype Manta DEL/DUP calls against CNVpytor single-cell data,
# summarise call distributions, and export filtered tables plus figures.
# Requirements: pandas, numpy, seaborn, matplotlib, cnvpytor; inputs in
# `samples_tissue`, `manta_vcfs/`, and `pytor_phased/`.

from __future__ import annotations

import gzip
import math
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

try:
    import cnvpytor
except ImportError:  # pragma: no cover - dependency may be absent during dry runs.
    cnvpytor = None


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SAMPLES_TISSUE = Path("samples_tissue") # file where each line is sample name space separated with tissue type (C - colon, L - lung)
VCF_DIR = Path("manta_vcfs")
PYTOR_DIR = Path("pytor")
SIZE_THRESH = 100_000
BIN_SIZE = 100_000
BLACKLIST = {} # samples to be removed
FIG_DIR = Path("figures")
OUTPUT_ALL = Path("manta_filtered_100k.csv")
OUTPUT_GENOTYPED = Path("manta_filtered_100k_genotyped.csv")

SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
END_RE = re.compile(r"(?:^|;)END=(\d+)")
SVLEN_RE = re.compile(r"(?:^|;)SVLEN=(-?\d+)")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_samples(path: Path) -> List[Tuple[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Required sample list {path} not found.")

    samples: List[Tuple[str, str]] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            sample = parts[0]
            tissue_code = parts[1].upper()
            samples.append((sample, tissue_code))
    return samples


def normalize_svtype(raw: Optional[str]) -> Optional[str]:
    if not raw:
        return None
    base = raw.split(":")[0].upper()
    if base.startswith("DEL"):
        return "DEL"
    if base.startswith("DUP"):
        return "DUP"
    if base.startswith("BND"):
        return "BND"
    return base


def parse_info(info: str) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    svtype = normalize_svtype(SVTYPE_RE.search(info).group(1)) if SVTYPE_RE.search(info) else None
    end = int(END_RE.search(info).group(1)) if END_RE.search(info) else None
    svlen = int(SVLEN_RE.search(info).group(1)) if SVLEN_RE.search(info) else None
    return svtype, end, svlen


def event_size(pos: int, end: Optional[int], svlen: Optional[int]) -> float:
    if svlen is not None and not math.isnan(svlen):
        return abs(int(svlen))
    if end is not None:
        return abs(int(end) - int(pos)) + 1
    return math.nan


def split_pair(val: str) -> Tuple[float, float]:
    try:
        ref, alt = val.split(",")
        return int(ref), int(alt)
    except Exception:
        return math.nan, math.nan


def parse_format_block(fmt: str, values: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    keys = fmt.split(":")
    vals = values.split(":")
    if len(vals) < len(keys):
        vals += [""] * (len(keys) - len(vals))
    for key, value in zip(keys, vals):
        key = key.strip()
        value = value.strip()
        if "," in value:
            ref, alt = split_pair(value)
            out[f"{key}_ref"] = ref
            out[f"{key}_alt"] = alt
        else:
            if value in {".", ""}:
                out[key] = math.nan
            else:
                try:
                    out[key] = int(value)
                except ValueError:
                    out[key] = value
    return out


def parse_vcf(sample: str, tissue: str, vcf_path: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    header_samples: Optional[Sequence[str]] = None

    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                header_samples = cols[9:11] if len(cols) >= 11 else None
                continue
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]
            fmt = cols[8]
            if filt != "PASS":
                continue

            svtype, end_v, svlen_v = parse_info(info)
            if svtype not in {"DEL", "DUP"}:
                continue

            pos_i = int(pos)
            size = event_size(pos_i, end_v, svlen_v)
            if not (math.isfinite(size) and size >= SIZE_THRESH):
                continue

            norm_vals = cols[9] if len(cols) >= 10 else ""
            sc_vals = cols[10] if len(cols) >= 11 else ""
            norm_dict = parse_format_block(fmt, norm_vals) if norm_vals else {}
            sc_dict = parse_format_block(fmt, sc_vals) if sc_vals else {}

            flattened = {f"{k}_norm": v for k, v in norm_dict.items()}
            flattened.update({f"{k}_sc": v for k, v in sc_dict.items()})

            row = {
                "cohort_sample": sample,
                "tissue": tissue,
                "vcf_normal_label": header_samples[0] if header_samples else "NORMAL_1",
                "vcf_sc_label": header_samples[1] if header_samples and len(header_samples) > 1 else "SC_2",
                "SVTYPE": svtype,
                "#CHROM": chrom,
                "POS": pos_i,
                "END": end_v,
                "SVLEN": svlen_v,
                "event_size": size,
                "ID": vid,
                "REF": ref,
                "ALT": alt,
                "QUAL": None if qual == "." else qual,
                "FILTER": filt,
                "INFO": info,
            }
            row.update(flattened)
            rows.append(row)
    return rows


def load_manta_calls(samples: Sequence[Tuple[str, str]]) -> Tuple[pd.DataFrame, List[str]]:
    rows: List[Dict[str, object]] = []
    missing: List[str] = []
    for sample, tissue in samples:
        vcf_path = VCF_DIR / f"{sample}.somaticSV.vcf.gz"
        if not vcf_path.exists():
            missing.append(sample)
            continue
        rows.extend(parse_vcf(sample, tissue, vcf_path))

    df = pd.DataFrame(rows)
    if df.empty:
        return df, missing

    meta_cols = [
        "cohort_sample",
        "tissue",
        "vcf_normal_label",
        "vcf_sc_label",
        "SVTYPE",
        "#CHROM",
        "POS",
        "END",
        "SVLEN",
        "event_size",
    ]
    count_cols_norm = sorted(col for col in df.columns if col.endswith("_norm"))
    count_cols_sc = sorted(col for col in df.columns if col.endswith("_sc"))
    vcf_cols = ["ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    ordered = meta_cols + count_cols_norm + count_cols_sc + vcf_cols
    df = df[[col for col in ordered if col in df.columns]]
    return df, missing


def summarise_counts(df: pd.DataFrame, samples: Sequence[Tuple[str, str]]) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["cohort_sample", "DEL", "DUP", "tissue", "Total_SVs"])

    all_samples = [sample for sample, _ in samples]
    counts = df.groupby(["cohort_sample", "SVTYPE"]).size().unstack(fill_value=0)
    counts = counts.reindex(index=all_samples, fill_value=0)
    counts = counts.reindex(columns=["DEL", "DUP"], fill_value=0)
    sample_to_tissue = dict(samples)
    counts["tissue"] = counts.index.map(sample_to_tissue.get)
    counts["Total_SVs"] = counts["DEL"] + counts["DUP"]
    return counts


def plot_histograms(counts: pd.DataFrame, output_path: Path) -> None:
    if counts.empty:
        return
    bins = np.arange(0, 22) - 0.5
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    axes[0].hist(counts["DEL"], bins=bins, edgecolor="black")
    axes[0].set_xticks(range(0, 21))
    axes[0].set_xlim(-0.5, 20.5)
    axes[0].set_title("Histogram of Deletions per Sample")
    axes[0].set_xlabel("# Deletions")
    axes[0].set_ylabel("# Samples")

    axes[1].hist(counts["DUP"], bins=bins, edgecolor="black")
    axes[1].set_xticks(range(0, 21))
    axes[1].set_xlim(-0.5, 20.5)
    axes[1].set_title("Histogram of Duplications per Sample")
    axes[1].set_xlabel("# Duplications")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Wrote histogram figure to {output_path}")


def plot_tissue_summary(df_filtered: pd.DataFrame, counts: pd.DataFrame) -> None:
    if df_filtered.empty or counts.empty:
        return

    agg = counts.groupby("tissue")[["DEL", "DUP"]].mean()
    fig, ax = plt.subplots(figsize=(6, 5))
    agg.plot(kind="bar", ax=ax)
    ax.set_title("Average #SVs per Sample by Tissue")
    ax.set_ylabel("Average count per sample")
    plt.xticks(rotation=0)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "tissue_avg_counts.png", dpi=200)
    plt.close(fig)
    print(f"Wrote tissue average count plot to {FIG_DIR / 'tissue_avg_counts.png'}")

    counts_plot = counts.copy()
    counts_plot["Total_SVs"] = counts_plot["DEL"] + counts_plot["DUP"]
    plt.figure(figsize=(6, 5))
    sns.violinplot(x="tissue", y="Total_SVs", data=counts_plot, inner="box", cut=0)
    plt.title("Distribution of SV Counts per Sample by Tissue")
    plt.ylabel("# SVs per sample")
    plt.xlabel("Tissue type")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "tissue_total_sv_violin.png", dpi=200)
    plt.close()
    print(f"Wrote tissue SV count violin to {FIG_DIR / 'tissue_total_sv_violin.png'}")

    plt.figure(figsize=(6, 5))
    sns.violinplot(x="tissue", y="event_size", data=df_filtered, inner="box", cut=0, scale="width")
    plt.title("Distribution of Event Sizes by Tissue")
    plt.ylabel("Event size (bp)")
    plt.xlabel("Tissue type")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "tissue_event_size_violin.png", dpi=200)
    plt.close()
    print(f"Wrote tissue event size violin to {FIG_DIR / 'tissue_event_size_violin.png'}")


class ViewerCache:
    def __init__(self, bin_size: int):
        self.bin_size = bin_size
        self.cache: Dict[str, cnvpytor.Viewer] = {}

    def get(self, sample: str) -> Optional[cnvpytor.Viewer]:
        if cnvpytor is None:
            return None
        if sample not in self.cache:
            pytor_path = PYTOR_DIR / f"{sample}.pytor"
            if not pytor_path.exists():
                return None
            self.cache[sample] = cnvpytor.Viewer([str(pytor_path)], {"bin_size": self.bin_size})
        return self.cache[sample]


def infer_end(pos: int, end: Optional[float], svlen: Optional[float]) -> float:
    if pd.notna(end):
        return float(end)
    if pd.notna(svlen):
        try:
            return pos + abs(float(svlen)) - 1
        except Exception:
            return math.nan
    return math.nan


def compute_rd(target_df: pd.DataFrame) -> pd.DataFrame:
    if target_df.empty:
        target_df["rd"] = pd.Series([], dtype=float)
        return target_df

    if cnvpytor is None:
        print("cnvpytor not available; skipping rd genotyping.")
        target_df["rd"] = math.nan
        return target_df

    cache = ViewerCache(BIN_SIZE)

    def row_rd(row: pd.Series) -> float:
        try:
            chrom = row["#CHROM"]
            pos = int(row["POS"])
            end = infer_end(pos, row.get("END", math.nan), row.get("SVLEN", math.nan))
            if not math.isfinite(end):
                return math.nan
            region = f"{chrom}:{pos}-{int(end)}"
            viewer = cache.get(str(row["cohort_sample"]))
            if viewer is None:
                return math.nan
            result = viewer.genotype([BIN_SIZE], region)
            return float(result[0][-1]) if result and result[0] else math.nan
        except Exception:
            return math.nan

    target_df["rd"] = target_df.apply(row_rd, axis=1)
    return target_df


def plot_rd_distributions(df: pd.DataFrame) -> None:
    if df.empty:
        return
    plot_df = df.dropna(subset=["rd"])
    if plot_df.empty:
        print("No rd values available for plotting.")
        return

    plt.figure(figsize=(6, 5))
    sns.violinplot(x="SVTYPE", y="rd", data=plot_df, inner="box", cut=0, scale="width")
    plt.title("Distribution of Copy Number (rd) by SV Type")
    plt.xlabel("SV Type")
    plt.ylabel("Copy number (rd)")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "rd_by_svtype_violin.png", dpi=200)
    plt.close()
    print(f"Wrote rd violin plot to {FIG_DIR / 'rd_by_svtype_violin.png'}")

    plt.figure(figsize=(7, 5))
    sns.kdeplot(data=plot_df, x="rd", hue="SVTYPE", common_norm=False, bw_adjust=0.5)
    plt.title("Density of Copy Number (rd) by SV Type")
    plt.xlabel("Copy number (rd)")
    plt.ylabel("Density")
    plt.xlim(0, 4)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "rd_by_svtype_kde.png", dpi=200)
    plt.close()
    print(f"Wrote rd density plot to {FIG_DIR / 'rd_by_svtype_kde.png'}")


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)

    samples = load_samples(SAMPLES_TISSUE)
    df, missing = load_manta_calls(samples)

    if missing:
        print(f"Missing VCF files for {len(missing)} samples: {', '.join(sorted(missing))}")

    if df.empty:
        print("No qualifying DEL/DUP PASS calls found. Exiting.")
        return

    print(f"Total DEL/DUP PASS calls (size ≥ {SIZE_THRESH:,}): {len(df)}")

    counts = summarise_counts(df, samples)
    plot_histograms(counts, FIG_DIR / "sv_counts_hist.png")

    counts_sorted = counts.sort_values("Total_SVs", ascending=False)
    print("SV counts per sample:")
    print(counts_sorted.to_string())

    df_filtered = df[~df["cohort_sample"].isin(BLACKLIST)].copy()
    print(f"Original df: {len(df)} rows")
    print(f"Filtered df: {len(df_filtered)} rows (removed {len(df) - len(df_filtered)} rows)")

    plot_tissue_summary(df_filtered, counts_sorted)

    target_df = df_filtered.copy()
    target_df = compute_rd(target_df)

    target_df.to_csv(OUTPUT_ALL, index=False)
    print(f"Wrote combined table with rd to {OUTPUT_ALL}")

    plot_rd_distributions(target_df)

    df_rd_filtered = target_df[(target_df["rd"] < 1.5) | (target_df["rd"] > 2.5)].copy()
    print(f"Original rows: {len(target_df)}")
    print(f"Filtered rows: {len(df_rd_filtered)}")
    df_rd_filtered.to_csv(OUTPUT_GENOTYPED, index=False)
    print(f"Wrote genotyped filtered table to {OUTPUT_GENOTYPED}")


if __name__ == "__main__":
    main()

# Script to aggregate CNVpytor call tables, classify large CNV events,
# generate genome-wide summary figures, and export filtered call statistics.
# Requirements: pandas, numpy, matplotlib; CNVpytor call tables in `calls/`.

from __future__ import annotations

import glob
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")  # Use a non-interactive backend for script execution.
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


SUFFIX = ".2d.100k.tsv"
PATH_PATTERN = Path("calls") / f"*{SUFFIX}"
COLUMN_NAMES = [
    "CNV_type",
    "CNV_region",
    "CNV_size",
    "CNV_level",
    "e-val1",
    "e-val2",
    "e-val3",
    "e-val4",
    "q0",
    "pN",
    "pNS",
    "pP",
    "bin_size",
    "n",
    "delta_BAF",
    "e-val1_rd_repeated",
    "baf_eval",
    "hets",
    "homs",
    "cn_1",
    "genotype_1",
    "likelihood_1",
    "cf_1",
    "cn_2",
    "genotype_2",
    "likelihood_2",
    "cf_2",
]

HG38_CHROM_SIZES = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chr3": 198_295_559,
    "chr4": 190_214_555,
    "chr5": 181_538_259,
    "chr6": 170_805_979,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr9": 138_394_717,
    "chr10": 133_797_422,
    "chr11": 135_086_622,
    "chr12": 133_275_309,
    "chr13": 114_364_328,
    "chr14": 107_043_718,
    "chr15": 101_991_189,
    "chr16": 90_338_345,
    "chr17": 83_257_441,
    "chr18": 80_373_285,
    "chr19": 58_617_616,
    "chr20": 64_444_167,
    "chr21": 46_709_983,
    "chr22": 50_818_468,
    "chrX": 156_040_895,
    "chrY": 57_227_415,
}

BIN_SIZE = 100_000
MAX_GAP_SIZE = 0  # Number of bins allowed when merging adjacent CNV segments.
SEGMENT_SIZE_THRESHOLD = 100  # Minimum number of bins to treat a segment as large.
FIG_DIR = Path("figures")


def load_cnv_tables(pattern: Path, suffix: str) -> pd.DataFrame:
    files = sorted(glob.glob(str(pattern)))
    if not files:
        print(f"No files found matching {pattern}.")
        return pd.DataFrame(columns=["sample", *COLUMN_NAMES])

    frames: List[pd.DataFrame] = []
    for file_path in files:
        sample_name = Path(file_path).name
        if sample_name.endswith(suffix):
            sample_name = sample_name[: -len(suffix)]
        df = pd.read_csv(file_path, sep="\t", header=None, names=COLUMN_NAMES)
        df.insert(0, "sample", sample_name)
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True)
    print(f"Loaded {len(files)} CNV tables with {combined.shape[0]} total rows.")
    return combined


def annotate_cnv_types(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    deletion_filter = (df["CNV_level"] < 0.75) & (df["n"] > 9) & (df["delta_BAF"].abs() > 0.4)
    duplication_filter = (
        (df["CNV_level"] > 1.3)
        & (df["n"] > 9)
        & (df["delta_BAF"].abs() > 0.08)
        & (df["delta_BAF"].abs() < 0.24)
    )
    cnnloh_filter = (
        (df["CNV_level"] > 0.8)
        & (df["CNV_level"] < 1.2)
        & (df["n"] > 9)
        & (df["delta_BAF"].abs() > 0.4)
    )

    conditions = [deletion_filter, duplication_filter, cnnloh_filter]
    labels = ["deletion", "duplication", "cnnloh"]
    df = df.copy()
    df["CNV_type2"] = np.select(conditions, labels, default=pd.NA)

    filtered = df[df["CNV_type2"].notna()].copy()
    print(f"Annotated CNV types. Found {len(filtered)} rows passing filters.")
    return df, filtered


def build_genome_index(
    chrom_sizes: Dict[str, int], bin_size: int
) -> Tuple[List[str], Dict[str, int], int]:
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    cumulative_bins: Dict[str, int] = {}
    total_bins = 0
    for chrom in chrom_order:
        cumulative_bins[chrom] = total_bins
        chrom_bins = int(math.ceil(chrom_sizes.get(chrom, 0) / bin_size))
        total_bins += chrom_bins
    return chrom_order, cumulative_bins, total_bins


def populate_plot_matrix(
    filtered_df: pd.DataFrame,
    chrom_order: List[str],
    cumulative_bins: Dict[str, int],
    chrom_sizes: Dict[str, int],
    bin_size: int,
    total_bins: int,
) -> pd.DataFrame:
    if filtered_df.empty:
        return pd.DataFrame()

    samples = sorted(filtered_df["sample"].unique())
    plot_matrix = pd.DataFrame(0, index=samples, columns=range(total_bins), dtype=int)
    cnv_to_int = {"deletion": 1, "duplication": 2, "cnnloh": 3}

    for _, row in filtered_df.iterrows():
        sample = row["sample"]
        cnv_type = row["CNV_type2"]
        region = row["CNV_region"]

        if pd.isna(cnv_type) or not isinstance(region, str) or ":" not in region:
            continue

        try:
            chrom, pos_range = region.split(":")
            start_pos, end_pos = map(int, pos_range.split("-"))
        except (ValueError, IndexError):
            continue

        if chrom not in cumulative_bins:
            continue

        start_bin = max((start_pos - 1) // bin_size, 0)
        end_bin = math.ceil(end_pos / bin_size)
        genome_offset = cumulative_bins[chrom]
        cnv_value = cnv_to_int.get(cnv_type)

        if cnv_value is None:
            continue

        chrom_bins = int(math.ceil(chrom_sizes.get(chrom, 0) / bin_size))
        for chrom_bin in range(start_bin, min(end_bin, chrom_bins)):
            genome_bin_index = genome_offset + chrom_bin
            if genome_bin_index < total_bins:
                plot_matrix.at[sample, genome_bin_index] = cnv_value

    return plot_matrix


def plot_heatmap(
    matrix: pd.DataFrame,
    chrom_order: List[str],
    cumulative_bins: Dict[str, int],
    output_path: Path,
    title: str,
) -> None:
    if matrix.empty:
        print(f"Skipping {title.lower()} – no data to plot.")
        return

    cmap = mcolors.ListedColormap(["#FFFFFF", "#E63946", "#0035ff", "#52B788"])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    fig_height = max(8, len(matrix.index) * 0.2)
    fig, ax = plt.subplots(figsize=(22, fig_height))
    ax.imshow(matrix.to_numpy(), aspect="auto", cmap=cmap, norm=norm, interpolation="none")

    ax.set_yticks(np.arange(len(matrix.index)))
    indexed_labels = [f"[{i}]: {sample}" for i, sample in enumerate(matrix.index)]
    ax.set_yticklabels(indexed_labels)
    ax.set_ylabel("Sample ID", fontweight="bold")

    xtick_locs = [cumulative_bins[c] for c in chrom_order]
    ax.set_xticks(xtick_locs)
    ax.set_xticklabels(chrom_order, rotation=90, ha="center", fontsize=9)
    ax.set_xlabel("Genomic Position (hg38)", fontweight="bold")

    for x in xtick_locs[1:]:
        ax.axvline(x - 0.5, color="lightgray", linestyle="--", linewidth=0.6)

    legend_elements = [
        Patch(facecolor="#E63946", edgecolor="black", label="Deletion"),
        Patch(facecolor="#0035ff", edgecolor="black", label="Duplication"),
        Patch(facecolor="#52B788", edgecolor="black", label="CNN-LOH"),
    ]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.01, 1), loc="upper left", title="CNV Type")

    ax.set_title(title, fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.95, 1])
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Wrote heatmap to {output_path}")


def merge_segments(
    matrix: pd.DataFrame,
    chrom_order: List[str],
    cumulative_bins: Dict[str, int],
    chrom_sizes: Dict[str, int],
    bin_size: int,
    total_bins: int,
    max_gap_size: int,
) -> pd.DataFrame:
    if matrix.empty:
        return matrix

    chromosome_map = np.empty(total_bins, dtype=object)
    for chrom in chrom_order:
        start = cumulative_bins[chrom]
        chrom_bins = int(math.ceil(chrom_sizes.get(chrom, 0) / bin_size))
        end = start + chrom_bins
        chromosome_map[start:end] = chrom

    merged = matrix.copy()
    for sample in merged.index:
        row = merged.loc[sample]
        for cnv_type in (1, 2, 3):
            indices = np.where(row.to_numpy() == cnv_type)[0]
            if len(indices) < 2:
                continue
            for left, right in zip(indices[:-1], indices[1:]):
                gap = right - left - 1
                same_chrom = chromosome_map[left] == chromosome_map[right]
                if 0 < gap <= max_gap_size and same_chrom:
                    merged.loc[sample, left:right] = cnv_type
    return merged


def calculate_stats(
    matrix: pd.DataFrame,
    chrom_order: List[str],
    cumulative_bins: Dict[str, int],
    chrom_sizes: Dict[str, int],
    bin_size: int,
    total_bins: int,
    segment_size_threshold: int,
) -> pd.DataFrame:
    stats: List[Dict[str, float]] = []
    for sample in matrix.index:
        row = matrix.loc[sample]
        cnv_blocks_all: List[pd.DataFrame] = []

        for chrom in chrom_order:
            start = cumulative_bins[chrom]
            chrom_bins = int(math.ceil(chrom_sizes.get(chrom, 0) / bin_size))
            end = start + chrom_bins
            chrom_row = row.iloc[start:end]
            if chrom_row.empty:
                continue

            blocks = (chrom_row != chrom_row.shift(1)).cumsum()
            block_stats = chrom_row.groupby(blocks).agg(["first", "size"])
            cnv_blocks = block_stats[block_stats["first"] != 0]
            if not cnv_blocks.empty:
                cnv_blocks_all.append(cnv_blocks)

        if cnv_blocks_all:
            cnv_blocks_sample = pd.concat(cnv_blocks_all)
        else:
            cnv_blocks_sample = pd.DataFrame(columns=["first", "size"])

        large_mask = cnv_blocks_sample["size"] > segment_size_threshold
        num_large_del = int(((cnv_blocks_sample["first"] == 1) & large_mask).sum())
        num_large_dup = int(((cnv_blocks_sample["first"] == 2) & large_mask).sum())
        num_large_loh = int(((cnv_blocks_sample["first"] == 3) & large_mask).sum())

        len_del = int(cnv_blocks_sample.loc[cnv_blocks_sample["first"] == 1, "size"].sum())
        len_dup = int(cnv_blocks_sample.loc[cnv_blocks_sample["first"] == 2, "size"].sum())
        len_loh = int(cnv_blocks_sample.loc[cnv_blocks_sample["first"] == 3, "size"].sum())

        total_large = num_large_del + num_large_dup + num_large_loh
        genome_altered = len_del + len_dup + len_loh
        percent_genome_altered = genome_altered / total_bins if total_bins else 0.0

        stats.append(
            {
                "sample": sample,
                "total_large_segments": total_large,
                "percent_genome_altered": percent_genome_altered,
                "large_deletions": num_large_del,
                "large_duplications": num_large_dup,
                "large_cnnloh": num_large_loh,
                "percent_del": len_del / total_bins if total_bins else 0.0,
                "percent_dup": len_dup / total_bins if total_bins else 0.0,
                "percent_loh": len_loh / total_bins if total_bins else 0.0,
            }
        )

    stats_df = pd.DataFrame(stats).set_index("sample")
    print("Calculated per-sample CNV statistics.")
    if not stats_df.empty:
        flagged = stats_df[
            (stats_df["total_large_segments"] > 10) | (stats_df["percent_genome_altered"] > 0.15)
        ]
        if not flagged.empty:
            print("Samples with substantial CNV burden:")
            print(flagged.to_string())
    return stats_df


def plot_summary_scatter(stats_df: pd.DataFrame, output_path: Path) -> None:
    if stats_df.empty:
        return

    fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
    ax.scatter(
        stats_df["percent_genome_altered"] * 100,
        stats_df["total_large_segments"],
    )
    ax.set_xlabel("Percentage of Genome Affected (%)")
    ax.set_ylabel("Number of Large Calls")
    ax.grid(True, linestyle="--", alpha=0.6)

    for idx, (sample, row) in enumerate(stats_df.iterrows()):
        ax.text(
            row["percent_genome_altered"] * 100,
            row["total_large_segments"],
            f"  [{idx}]",
            fontsize=5,
        )

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Wrote scatter plot to {output_path}")


def plot_summary_bar(stats_df: pd.DataFrame, output_path: Path) -> None:
    if stats_df.empty:
        return

    totals = {
        "Deletion": int(stats_df["large_deletions"].sum()),
        "Duplication": int(stats_df["large_duplications"].sum()),
        "CNN-LOH": int(stats_df["large_cnnloh"].sum()),
    }
    colors = ["#E63946", "#0035ff", "#52B788"]

    fig, ax = plt.subplots(figsize=(7, 5), dpi=100)
    bars = ax.bar(totals.keys(), totals.values(), color=colors, edgecolor="black")
    ax.set_xlabel("CNV Type", fontweight="bold")
    ax.set_ylabel("Total Number of Large Segments", fontweight="bold")
    ax.set_title(
        f"Total Count of Large CNV Segments by Type (>{SEGMENT_SIZE_THRESHOLD} bins)",
        fontweight="bold",
    )
    ax.grid(axis="y", linestyle="--", alpha=0.7)
    ax.set_ylim(0, max(totals.values()) * 1.15 if totals else 1)

    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2.0, height, int(height), va="bottom", ha="center")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Wrote bar plot to {output_path}")


def read_tissue_map(path: Path) -> Dict[str, str]:
    if not path.exists():
        print(f"Tissue map {path} not found; skipping tissue-stratified plots.")
        return {}

    mapping: Dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("---") or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            tissue_code = parts[1].upper()
            tissue = "Colon" if tissue_code.startswith("C") else ("Lung" if tissue_code.startswith("L") else tissue_code)
            mapping[parts[0]] = tissue
    return mapping


def plot_tissue_counts(
    stats_df: pd.DataFrame,
    tissue_map: Dict[str, str],
    output_prefix: Path,
) -> None:
    if stats_df.empty or not tissue_map:
        return

    columns_needed = [
        "large_deletions",
        "large_duplications",
        "large_cnnloh",
        "total_large_segments",
        "percent_genome_altered",
        "percent_del",
        "percent_dup",
        "percent_loh",
    ]
    missing_cols = [col for col in columns_needed if col not in stats_df.columns]
    if missing_cols:
        print(f"Missing columns {missing_cols}; skipping tissue counts.")
        return

    full_samples = sorted(tissue_map)
    stats_t = stats_df.reindex(full_samples).fillna(0)
    stats_t.index.name = "sample"
    stats_t["tissue"] = stats_t.index.map(tissue_map.get)
    stats_t = stats_t[stats_t["tissue"].isin(["Colon", "Lung"])]

    for col in ["large_deletions", "large_duplications", "large_cnnloh", "total_large_segments"]:
        stats_t[col] = stats_t[col].astype(int)

    if stats_t.empty:
        print("No overlapping samples between stats and tissue map; skipping tissue plots.")
        return

    cnv_types = ["Deletions", "Duplications", "CNN-LOH"]
    cnv_cols = {
        "Deletions": "large_deletions",
        "Duplications": "large_duplications",
        "CNN-LOH": "large_cnnloh",
    }
    cnv_colors = {"Deletions": "#E63946", "Duplications": "#0035ff", "CNN-LOH": "#52B788"}
    tissues = ["Colon", "Lung"]

    counts_by_tissue_type = []
    for tissue in tissues:
        df_t = stats_t[stats_t["tissue"] == tissue]
        counts_by_tissue_type.append(
            [int((df_t[cnv_cols[name]] > 0).sum()) for name in cnv_types]
        )

    n_types = len(cnv_types)
    group_gap = 1
    group_width = n_types + group_gap
    x_positions: List[int] = []
    colors: List[str] = []
    for idx, tissue in enumerate(tissues):
        base = idx * group_width
        for offset, cnv_type in enumerate(cnv_types):
            x_positions.append(base + offset)
            colors.append(cnv_colors[cnv_type])

    values = [counts_by_tissue_type[g][i] for g in range(len(tissues)) for i in range(n_types)]

    fig, ax = plt.subplots(figsize=(10, 5), dpi=100)
    bars = ax.bar(x_positions, values, edgecolor="black", color=colors)
    ax.set_xticks(x_positions, [cnv_types[i % n_types] for i in range(len(x_positions))])
    ax.set_ylabel("Number of Samples", fontweight="bold")
    ax.set_title(
        f"Samples with ≥1 Large CNV per Type, grouped by Tissue (>{SEGMENT_SIZE_THRESHOLD} bins)",
        fontweight="bold",
    )
    ax.grid(axis="y", linestyle="--", alpha=0.7)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim(0, max(values + [1]) * 1.15)

    for g, tissue in enumerate(tissues):
        center = g * group_width + (n_types - 1) / 2
        ax.text(
            center,
            -0.12,
            tissue,
            ha="center",
            va="top",
            transform=ax.get_xaxis_transform(),
            fontsize=10,
            fontweight="bold",
        )

    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height, f"{height}", va="bottom", ha="center", fontweight="bold")

    fig.subplots_adjust(bottom=0.22)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix("_cnv_types.png"), dpi=200)
    plt.close(fig)
    print(f"Wrote tissue CNV type plot to {output_prefix.with_suffix('_cnv_types.png')}")

    any_labels = ["Any Large CNV", "No Large CNV"]
    any_colors = {"Any Large CNV": "#1B9E77", "No Large CNV": "#CCCCCC"}

    counts_any = []
    for tissue in tissues:
        df_t = stats_t[stats_t["tissue"] == tissue]
        any_large = int((df_t["total_large_segments"] > 0).sum())
        none_large = int((df_t["total_large_segments"] == 0).sum())
        counts_any.append([any_large, none_large])

    n_types2 = len(any_labels)
    group_width2 = n_types2 + group_gap
    x2_positions: List[int] = []
    colors2: List[str] = []
    for idx, tissue in enumerate(tissues):
        base = idx * group_width2
        for offset, label in enumerate(any_labels):
            x2_positions.append(base + offset)
            colors2.append(any_colors[label])

    values2 = [counts_any[g][i] for g in range(len(tissues)) for i in range(n_types2)]

    fig, ax = plt.subplots(figsize=(9, 5), dpi=100)
    bars2 = ax.bar(x2_positions, values2, edgecolor="black", color=colors2)
    ax.set_xticks(x2_positions, [any_labels[i % n_types2] for i in range(len(x2_positions))])
    ax.set_ylabel("Number of Samples", fontweight="bold")
    ax.set_title(
        f"Any Large CNV vs None, grouped by Tissue (>{SEGMENT_SIZE_THRESHOLD} bins)",
        fontweight="bold",
    )
    ax.grid(axis="y", linestyle="--", alpha=0.7)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim(0, max(values2 + [1]) * 1.15)

    for g, tissue in enumerate(tissues):
        center = g * group_width2 + (n_types2 - 1) / 2
        ax.text(
            center,
            -0.12,
            tissue,
            ha="center",
            va="top",
            transform=ax.get_xaxis_transform(),
            fontsize=10,
            fontweight="bold",
        )

    for bar in bars2:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height, f"{height}", va="bottom", ha="center", fontweight="bold")

    fig.subplots_adjust(bottom=0.22)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix("_any_cnv.png"), dpi=200)
    plt.close(fig)
    print(f"Wrote tissue any/none plot to {output_prefix.with_suffix('_any_cnv.png')}")


def save_large_calls(
    matrix: pd.DataFrame,
    chrom_order: List[str],
    cumulative_bins: Dict[str, int],
    chrom_sizes: Dict[str, int],
    bin_size: int,
    tissue_map: Dict[str, str],
    segment_size_threshold: int,
    output_path: Path,
) -> None:
    if matrix.empty:
        print("No merged CNV matrix available; skipping large call export.")
        return

    int_to_cnv_type = {1: "deletion", 2: "duplication", 3: "cnnloh"}
    records: List[Dict[str, object]] = []

    for sample in matrix.index:
        row = matrix.loc[sample]
        for chrom in chrom_order:
            start_abs = cumulative_bins[chrom]
            chrom_bins = int(math.ceil(chrom_sizes.get(chrom, 0) / bin_size))
            end_abs = start_abs + chrom_bins
            chrom_row = row.iloc[start_abs:end_abs]
            if chrom_row.empty:
                continue

            blocks = (chrom_row != chrom_row.shift(1)).cumsum()
            for block_id, block_series in chrom_row.groupby(blocks):
                cnv_type_int = int(block_series.iloc[0])
                if cnv_type_int == 0:
                    continue

                block_len = int(block_series.size)
                if block_len <= segment_size_threshold:
                    continue

                rel_indices = block_series.index - start_abs
                start_pos = int(rel_indices.min()) * bin_size + 1
                end_pos = int(rel_indices.max() + 1) * bin_size
                cnv_size = end_pos - start_pos + 1

                records.append(
                    {
                        "sample": sample,
                        "tissue_type": tissue_map.get(sample, "NA"),
                        "region": f"{chrom}:{start_pos}-{end_pos}",
                        "cnv_type": int_to_cnv_type.get(cnv_type_int, f"type_{cnv_type_int}"),
                        "cnv_size": cnv_size,
                        "size_of_chromosome": int(chrom_sizes.get(chrom, 0)),
                    }
                )

    if not records:
        print("No large CNV calls identified; skipping export.")
        return

    large_calls_df = (
        pd.DataFrame.from_records(records)
        .sort_values(["sample", "region", "cnv_type"], ignore_index=True)
    )
    large_calls_df.to_csv(output_path, sep="\t", index=False)
    print(f"Saved {len(large_calls_df)} large calls to {output_path}")


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)

    combined_df = load_cnv_tables(PATH_PATTERN, SUFFIX)
    if combined_df.empty:
        return

    annotated_df, filtered_df = annotate_cnv_types(combined_df)
    print(f"Combined table shape: {annotated_df.shape}")
    print(f"Filtered table shape: {filtered_df.shape}")

    chrom_order, cumulative_bins, total_bins = build_genome_index(HG38_CHROM_SIZES, BIN_SIZE)
    plot_matrix = populate_plot_matrix(
        filtered_df, chrom_order, cumulative_bins, HG38_CHROM_SIZES, BIN_SIZE, total_bins
    )

    plot_heatmap(
        plot_matrix,
        chrom_order,
        cumulative_bins,
        FIG_DIR / "genomewide_cnv_calls.png",
        "Genome-wide CNV Calls by Sample",
    )

    merged_matrix = merge_segments(
        plot_matrix,
        chrom_order,
        cumulative_bins,
        HG38_CHROM_SIZES,
        BIN_SIZE,
        total_bins,
        MAX_GAP_SIZE,
    )
    plot_heatmap(
        merged_matrix,
        chrom_order,
        cumulative_bins,
        FIG_DIR / "genomewide_cnv_calls_merged.png",
        "Genome-wide CNV Calls (Segments Merged)",
    )

    stats_df = calculate_stats(
        merged_matrix, chrom_order, cumulative_bins, HG38_CHROM_SIZES, BIN_SIZE, total_bins, SEGMENT_SIZE_THRESHOLD
    )
    if stats_df.empty:
        return

    plot_summary_scatter(stats_df, FIG_DIR / "summary_scatter.png")
    plot_summary_bar(stats_df, FIG_DIR / "summary_bar.png")

    tissue_map = read_tissue_map(Path("samples_tissue"))
    plot_tissue_counts(stats_df, tissue_map, FIG_DIR / "tissue_summary")

    save_large_calls(
        merged_matrix,
        chrom_order,
        cumulative_bins,
        HG38_CHROM_SIZES,
        BIN_SIZE,
        tissue_map,
        SEGMENT_SIZE_THRESHOLD,
        Path("large_calls.tsv"),
    )


if __name__ == "__main__":
    main()

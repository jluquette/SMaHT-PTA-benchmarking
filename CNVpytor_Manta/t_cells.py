#!/usr/bin/env python3
# Script to aggregate Manta SV calls at T-cell receptor loci, extract CNVpytor
# read-depth profiles, fit discrete copy-number states, and render summary
# figures. Requirements: pandas, numpy, matplotlib, cnvpytor; input files in
# `samples_tissue`, `manta_vcfs/`, and `pytor_phased/`.

from __future__ import annotations

import gzip
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

try:
    import cnvpytor  # type: ignore
except ImportError:  # pragma: no cover - dependency may be missing in some environments.
    cnvpytor = None


SAMPLES_TISSUE = Path("../samples_tissue")
VCF_DIR = Path("../manta_vcfs")
PYTOR_DIR = Path("../pytor")
FIG_DIR = Path("figures")
MIN_SV_SIZE = 10_000
PLOT_EVENT_SIZE_MIN = 20_000
# Update `CONTROL_SAMPLE` and `T_CELLS` with real sample identifiers before running.
BIN_SIZE = 10_000
CONTROL_SAMPLE = "control_sample"
T_CELLS = [
    "t_cell_sample_1",
    "t_cell_sample_2",
    "t_cell_sample_3",
    "t_cell_sample_4",
    "t_cell_sample_5",
]
T_CELL_LOCI = ["chr14:21621904-22552132", "chr7:142299011-142813287", "chr7:38240024-38368055"]
CN_MAX = 3
MAX_SWITCHES = 5


def parse_info(info_str: str) -> Dict[str, str]:
    info: Dict[str, str] = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
    return info


def parse_alt_bnd(alt_str: str) -> Tuple[Optional[str], Optional[int]]:
    match = re.search(r"([\[\]])(chr[\w\d]+):(\d+)([\[\]])", alt_str)
    if not match:
        return None, None
    _, chrom, pos, _ = match.groups()
    return chrom, int(pos)


def parse_sample_format(format_str: str, sample_data_str: str) -> Dict[str, float]:
    reads = {"PR_ref": 0, "PR_alt": 0, "SR_ref": 0, "SR_alt": 0}
    keys = format_str.split(":")
    values = sample_data_str.split(":")
    data_map = dict(zip(keys, values))

    if "PR" in data_map:
        try:
            ref, alt = map(int, data_map["PR"].split(","))
            reads["PR_ref"] = ref
            reads["PR_alt"] = alt
        except (ValueError, IndexError):
            pass
    if "SR" in data_map:
        try:
            ref, alt = map(int, data_map["SR"].split(","))
            reads["SR_ref"] = ref
            reads["SR_alt"] = alt
        except (ValueError, IndexError):
            pass
    return reads


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
            samples.append((parts[0], parts[1].upper()))
    return samples


def parse_vcf(sample: str, tissue: str, vcf_path: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="ignore") as fh:
        header_samples: Optional[Sequence[str]] = None
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

            chrom, pos, vid, ref, alt, qual, filt, info_str = cols[:8]
            fmt = cols[8]
            if filt != "PASS":
                continue

            info = parse_info(info_str)
            svtype = info.get("SVTYPE")
            if svtype not in {"DEL", "DUP", "BND"}:
                continue

            pos_i = int(pos)
            end_v = int(info["END"]) if "END" in info else None
            svlen_raw = info.get("SVLEN")
            svlen_v = int(svlen_raw) if svlen_raw is not None else None

            if svtype in {"DEL", "DUP"}:
                size = abs(svlen_v) if svlen_v is not None else (abs(int(end_v) - pos_i) if end_v else None)
                if size is None or size <= MIN_SV_SIZE:
                    continue
                end_chrom = chrom
                end_pos = end_v if end_v is not None else pos_i + size
            else:
                mate_chrom, mate_pos = parse_alt_bnd(alt)
                if mate_chrom is None:
                    continue
                if chrom == mate_chrom and abs(pos_i - mate_pos) <= MIN_SV_SIZE:
                    continue
                end_chrom, end_pos = mate_chrom, mate_pos
                size = abs(end_pos - pos_i)

            norm_reads = parse_sample_format(fmt, cols[9]) if len(cols) >= 10 else {}
            sc_reads = parse_sample_format(fmt, cols[10]) if len(cols) >= 11 else {}

            row = {
                "sample": sample,
                "tissue": tissue,
                "SVTYPE": svtype,
                "chrom": chrom,
                "pos": pos_i,
                "end_chrom": end_chrom,
                "end_pos": end_pos,
                "SVLEN": svlen_v,
                "event_size": size if size is not None else math.nan,
                "ID": vid,
                "info": info_str,
                "alt": alt,
                "qual": None if qual == "." else qual,
                "filter": filt,
                "total_alt_reads": sc_reads.get("PR_alt", 0) + sc_reads.get("SR_alt", 0),
                "normal_PR_ref": norm_reads.get("PR_ref"),
                "normal_PR_alt": norm_reads.get("PR_alt"),
                "normal_SR_ref": norm_reads.get("SR_ref"),
                "normal_SR_alt": norm_reads.get("SR_alt"),
                "sc_PR_ref": sc_reads.get("PR_ref"),
                "sc_PR_alt": sc_reads.get("PR_alt"),
                "sc_SR_ref": sc_reads.get("SR_ref"),
                "sc_SR_alt": sc_reads.get("SR_alt"),
            }
            rows.append(row)
    return rows


def load_sv_dataframe(samples: Sequence[Tuple[str, str]]) -> Tuple[pd.DataFrame, List[str]]:
    rows: List[Dict[str, object]] = []
    missing: List[str] = []
    for sample, tissue in samples:
        vcf_path = VCF_DIR / f"{sample}.somaticSV.vcf.gz"
        if not vcf_path.exists():
            missing.append(sample)
            continue
        rows.extend(parse_vcf(sample, tissue, vcf_path))
    return pd.DataFrame(rows), missing


class ViewerCache:
    def __init__(self, bin_size: int):
        self.bin_size = bin_size
        self._cache: Dict[str, cnvpytor.Viewer] = {}

    def get(self, sample: str) -> Optional[cnvpytor.Viewer]:
        if cnvpytor is None:
            return None
        if sample not in self._cache:
            pytor_path = PYTOR_DIR / f"{sample}.pytor"
            if not pytor_path.exists():
                return None
            self._cache[sample] = cnvpytor.Viewer([str(pytor_path)], {"bin_size": self.bin_size})
        return self._cache[sample]


def get_rd(sample: str, region: str, cache: ViewerCache, bin_size: int = BIN_SIZE) -> np.ndarray:
    viewer = cache.get(sample)
    if viewer is None:
        return np.array([])

    chrom, coords = region.split(":")
    start_str, end_str = coords.split("-")
    start, end = int(start_str), int(end_str)
    start_bin = start // bin_size
    end_bin = end // bin_size

    try:
        mean, _ = viewer.io[0].rd_normal_level(bin_size, cnvpytor.FLAG_GC_CORR)
        signal = viewer.io[0].get_signal(chrom, bin_size, "RD", cnvpytor.FLAG_GC_CORR)
    except Exception:
        return np.array([])

    segment = np.array(signal[start_bin : end_bin + 1], dtype=float)
    if segment.size == 0 or mean in (0, None):
        return np.array([])
    return segment / mean * 2.0


def solve_fitting_problem(read_depths: Sequence[float], max_switches: int = MAX_SWITCHES, max_cn: int = CN_MAX) -> List[int]:
    signal = np.asarray(read_depths, dtype=float)
    n = signal.size
    if n == 0:
        return []

    cost = np.full((n, n), np.inf)
    cost_val = np.zeros((n, n), dtype=int)

    prefix_sum = np.concatenate(([0.0], np.cumsum(signal)))
    prefix_sum_sq = np.concatenate(([0.0], np.cumsum(signal**2)))

    for i in range(n):
        for j in range(i, n):
            seg_len = j - i + 1
            seg_sum = prefix_sum[j + 1] - prefix_sum[i]
            seg_sq_sum = prefix_sum_sq[j + 1] - prefix_sum_sq[i]

            best_err = np.inf
            best_val = 0
            for val in range(max_cn + 1):
                err = seg_sq_sum - 2 * val * seg_sum + seg_len * (val**2)
                if err < best_err:
                    best_err = err
                    best_val = val
            cost[i][j] = best_err
            cost_val[i][j] = best_val

    dp = np.full((max_switches + 1, n + 1), np.inf)
    path = np.zeros((max_switches + 1, n + 1), dtype=int)
    dp[:, 0] = 0.0
    for i in range(1, n + 1):
        dp[0][i] = cost[0][i - 1]

    for k in range(1, max_switches + 1):
        for i in range(1, n + 1):
            best_cost = np.inf
            best_j = 0
            for j in range(1, i + 1):
                current = dp[k - 1][j - 1] + cost[j - 1][i - 1]
                if current < best_cost:
                    best_cost = current
                    best_j = j - 1
            dp[k][i] = best_cost
            path[k][i] = best_j

    final_errors = dp[:, n]
    best_k = int(np.argmin(final_errors))

    output = np.zeros(n, dtype=int)
    end_idx = n
    k = best_k
    while end_idx > 0:
        start_idx = path[k][end_idx]
        output[start_idx:end_idx] = cost_val[start_idx][end_idx - 1]
        end_idx = start_idx
        if k > 0:
            k -= 1
    return output.tolist()


def compute_rd_tracks(samples: Sequence[str], loci: Sequence[str], cache: ViewerCache) -> Dict[Tuple[str, str], Dict[str, np.ndarray]]:
    rd_data: Dict[Tuple[str, str], Dict[str, np.ndarray]] = {}
    for sample in samples:
        for locus in loci:
            signal = get_rd(sample, locus, cache, BIN_SIZE)
            fit = np.array(solve_fitting_problem(signal, MAX_SWITCHES, CN_MAX)) if signal.size else np.array([])
            rd_data[(sample, locus)] = {"raw": signal, "fit": fit}
    return rd_data


def filter_svs_for_region(
    sv_df: pd.DataFrame,
    sample: str,
    chrom: str,
    start: int,
    end: int,
    min_size: int = PLOT_EVENT_SIZE_MIN,
) -> pd.DataFrame:
    if sv_df.empty:
        return sv_df
    mask = (
        (sv_df["sample"] == sample)
        & (sv_df["chrom"] == chrom)
        & (sv_df["end_chrom"] == chrom)
        & (sv_df["pos"] >= start)
        & (sv_df["end_pos"] <= end)
        & ((sv_df["end_pos"] - sv_df["pos"]) > min_size)
    )
    return sv_df.loc[mask].copy()


def plot_step_grid(
    sv_df: pd.DataFrame,
    rd_tracks: Dict[Tuple[str, str], Dict[str, np.ndarray]],
    samples: Sequence[str],
    loci: Sequence[str],
    output_path: Path,
) -> None:
    if not samples or not loci:
        return
    fig, axes = plt.subplots(
        nrows=len(samples),
        ncols=len(loci),
        figsize=(3.0 * len(loci), 1.2 * len(samples)),
        squeeze=False,
    )

    for i, sample in enumerate(samples):
        for j, locus in enumerate(loci):
            ax = axes[i][j]
            match = re.match(r"(chr\w+):(\d+)-(\d+)", locus)
            if not match:
                continue
            chrom, start_str, end_str = match.groups()
            start, end = int(start_str), int(end_str)
            region_len = end - start
            padding = region_len // 2
            zoom_start = max(1, start - padding)
            zoom_end = end + padding
            zoom_locus = f"{chrom}:{zoom_start}-{zoom_end}"

            track = rd_tracks.get((sample, zoom_locus))
            if track is None:
                track = rd_tracks.get((sample, locus))
            if track is None:
                track = {"raw": np.array([]), "fit": np.array([])}

            raw = track["raw"]
            fit = track["fit"]
            if raw.size:
                ax.step(range(raw.size), raw, color="grey", linewidth=0.8)
            if fit.size:
                ax.step(range(fit.size), fit, color="black", linewidth=1.0)

            start_bin = (start - zoom_start) / BIN_SIZE
            end_bin = (end - zoom_start) / BIN_SIZE
            ax.axvline(start_bin, color="#2ca02c", linestyle="-", linewidth=1.2)
            ax.axvline(end_bin, color="#2ca02c", linestyle="-", linewidth=1.2)

            sv_subset = filter_svs_for_region(
                sv_df, sample, chrom, zoom_start, zoom_end, PLOT_EVENT_SIZE_MIN
            )
            if not sv_subset.empty:
                for idx, (_, row) in enumerate(sv_subset.iterrows()):
                    x_start = (row["pos"] - zoom_start) / BIN_SIZE
                    x_end = (row["end_pos"] - zoom_start) / BIN_SIZE
                    y_pos = 3.5 - idx * 0.3
                    ax.annotate(
                        "",
                        xy=(x_end, y_pos),
                        xytext=(x_start, y_pos),
                        arrowprops=dict(arrowstyle="<->", color="red", lw=1.0),
                    )

            ax.set_ylim(0, 4)
            ax.set_xticks([])
            ax.set_yticks([0, 1, 2, 3, 4])
            ax.set_yticklabels([])
            ax.grid(axis="y", linestyle="dotted", color="lightgray")
            ax.tick_params(axis="y", length=0)
            for spine in ax.spines.values():
                spine.set_visible(False)

    fig.subplots_adjust(left=0.08, right=0.92, top=0.96, bottom=0.08, wspace=0.04, hspace=0.08)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote step plot grid to {output_path}")


def build_cmap() -> LinearSegmentedColormap:
    vmax = 4.0
    nodes = [
        (0.0 / vmax, "white"),
        (0.1 / vmax, "white"),
        (0.9 / vmax, "#777777"),
        (1.1 / vmax, "#777777"),
        (1.9 / vmax, "#000000"),
        (2.1 / vmax, "#000000"),
        (2.9 / vmax, "#000077"),
        (3.1 / vmax, "#000077"),
        (4.0 / vmax, "#007700"),
    ]
    return LinearSegmentedColormap.from_list("custom_rd_cmap", nodes)


def plot_heatmap_grid(
    sv_df: pd.DataFrame,
    rd_tracks: Dict[Tuple[str, str], Dict[str, np.ndarray]],
    samples: Sequence[str],
    loci: Sequence[str],
    output_path: Path,
    use_fit: bool = True,
) -> None:
    if not samples or not loci:
        return
    cmap = build_cmap()

    heights: List[float] = []
    for _ in samples:
        heights.extend([4, 1.5])

    fig = plt.figure(figsize=(3.5 * len(loci), 1.3 * len(samples)))
    gs = gridspec.GridSpec(
        nrows=2 * len(samples),
        ncols=len(loci),
        figure=fig,
        wspace=0.1,
        hspace=0.15,
        height_ratios=heights,
    )

    im = None
    for i, sample in enumerate(samples):
        for j, locus in enumerate(loci):
            match = re.match(r"(chr\w+):(\d+)-(\d+)", locus)
            if not match:
                continue
            chrom, start_str, end_str = match.groups()
            start, end = int(start_str), int(end_str)
            region_len = end - start
            padding = region_len // 2
            zoom_start = max(1, start - padding)
            zoom_end = end + padding
            zoom_locus = f"{chrom}:{zoom_start}-{zoom_end}"

            track = rd_tracks.get((sample, zoom_locus))
            if track is None:
                track = rd_tracks.get((sample, locus))
            if track is None:
                track = {"raw": np.array([]), "fit": np.array([])}

            data = track["fit"] if use_fit else track["raw"]
            ax_map = fig.add_subplot(gs[2 * i, j])
            ax_sv = fig.add_subplot(gs[2 * i + 1, j], sharex=ax_map)

            if data.size == 0:
                data = np.zeros(1)
            im = ax_map.imshow(
                data.reshape(1, -1),
                cmap=cmap,
                aspect="auto",
                interpolation="nearest",
                vmin=0,
                vmax=4,
            )

            start_bin = (start - zoom_start) / BIN_SIZE
            end_bin = (end - zoom_start) / BIN_SIZE
            ax_map.axvline(start_bin, color="#da7400", linestyle=":", linewidth=2)
            ax_map.axvline(end_bin, color="#da7400", linestyle=":", linewidth=2)
            ax_map.set_yticks([])
            plt.setp(ax_map.get_xticklabels(), visible=False)
            ax_map.tick_params(axis="x", length=0)
            for spine in ax_map.spines.values():
                spine.set_visible(False)

            sv_subset = filter_svs_for_region(
                sv_df, sample, chrom, zoom_start, zoom_end, PLOT_EVENT_SIZE_MIN
            )
            ax_sv.set_yticks([])
            ax_sv.set_xticks([])
            ax_sv.set_ylim(0, max(len(sv_subset), 1))
            for spine in ax_sv.spines.values():
                spine.set_visible(False)

            if not sv_subset.empty:
                for idx, (_, row) in enumerate(sv_subset.iterrows()):
                    x_start = (row["pos"] - zoom_start) / BIN_SIZE
                    x_end = (row["end_pos"] - zoom_start) / BIN_SIZE
                    ax_sv.annotate(
                        "",
                        xy=(x_end, idx + 0.5),
                        xytext=(x_start, idx + 0.5),
                        arrowprops=dict(arrowstyle="<->", color="red", lw=1.3),
                    )

    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.7])
    if im is not None:
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label("Read Depth", rotation=-90, va="bottom")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote heatmap grid to {output_path}")


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    samples = load_samples(SAMPLES_TISSUE)
    sv_df, missing_vcfs = load_sv_dataframe(samples)
    if missing_vcfs:
        print(f"Missing VCF files for {len(missing_vcfs)} samples: {', '.join(sorted(missing_vcfs))}")

    summary_samples = [CONTROL_SAMPLE] + T_CELLS
    present = {sample for sample, _ in samples}
    absent = [s for s in summary_samples if s not in present]
    if absent:
        print(f"Warning: samples not found in tissue list: {', '.join(absent)}")

    plot_samples = [s for s in summary_samples if s in present]
    if cnvpytor is None:
        print("cnvpytor not available; plots will omit read-depth data.")

    cache = ViewerCache(BIN_SIZE)
    loci_zoom: List[str] = []
    for locus in T_CELL_LOCI:
        match = re.match(r"(chr\w+):(\d+)-(\d+)", locus)
        if match:
            chrom, start_str, end_str = match.groups()
            start, end = int(start_str), int(end_str)
            pad = (end - start) // 2
            zoom_start = max(1, start - pad)
            zoom_end = end + pad
            loci_zoom.append(f"{chrom}:{zoom_start}-{zoom_end}")
    combined_loci = list(dict.fromkeys(loci_zoom + T_CELL_LOCI))
    rd_tracks = compute_rd_tracks(plot_samples, combined_loci, cache)

    if sv_df.empty:
        print("No SV records found; figures will still be generated from RD tracks if available.")

    step_plot_path = FIG_DIR / "figC_step.png"
    plot_step_grid(sv_df, rd_tracks, plot_samples, T_CELL_LOCI, step_plot_path)

    heatmap_fit_path = FIG_DIR / "figC_heatmap_fit.png"
    plot_heatmap_grid(sv_df, rd_tracks, plot_samples, T_CELL_LOCI, heatmap_fit_path, use_fit=True)

    heatmap_raw_path = FIG_DIR / "figC_heatmap_raw.png"
    plot_heatmap_grid(sv_df, rd_tracks, plot_samples, T_CELL_LOCI, heatmap_raw_path, use_fit=False)


if __name__ == "__main__":
    main()

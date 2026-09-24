#!/usr/bin/env python3
"""Render Figure 3 from the validated raw 250-kb Cooler matrix.

Values are raw observed Hi-C contacts transformed as log10(raw count + 1).
Display limits come from the supplied parameter TSV. Cooler bin start/end
coordinates are converted to true cumulative sequence coordinates, so terminal
partial bins are displayed at their true widths.
"""
from __future__ import annotations
import argparse, csv, hashlib
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from common import set_pub_style, save_triplet_fixed


def read_parameters(path: Path):
    out = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames != ["parameter", "value"]:
            raise SystemExit(f"Expected parameter/value TSV: {path}")
        for row in reader:
            out[row["parameter"]] = row["value"]
    required = ["resolution_bp", "normalization", "transform", "vmin_log10",
                "vmax_log10", "coordinate_system", "primary_total_bp"]
    missing = [x for x in required if x not in out]
    if missing:
        raise SystemExit(f"Missing display parameters: {missing}")
    return out


def sha256_file(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def coord_formatter(total_mb):
    def fmt(x, _pos):
        if abs(x-total_mb) < 0.15:
            return f"{total_mb:.1f}"
        if abs(x-round(x)) < 1e-6:
            return str(int(round(x)))
        return f"{x:g}"
    return FuncFormatter(fmt)


def true_cumulative_edges(c):
    bins = c.bins()[:][["chrom", "start", "end"]].copy()
    chromnames = list(c.chromnames)
    chromsizes = c.chromsizes
    offsets, cum = {}, 0
    for chrom in chromnames:
        offsets[chrom] = cum
        cum += int(chromsizes[chrom])
    starts = bins["start"].to_numpy(dtype=np.int64) + bins["chrom"].map(offsets).to_numpy(dtype=np.int64)
    ends = bins["end"].to_numpy(dtype=np.int64) + bins["chrom"].map(offsets).to_numpy(dtype=np.int64)
    if len(starts) == 0 or starts[0] != 0:
        raise SystemExit("Unexpected or empty Cooler bin table")
    if not np.array_equal(starts[1:], ends[:-1]):
        bad = np.flatnonzero(starts[1:] != ends[:-1])[:10]
        raise SystemExit(f"Non-contiguous cumulative bin edges at indices {bad.tolist()}")
    edges = np.concatenate(([starts[0]], ends)).astype(np.float64)
    if int(edges[-1]) != int(cum):
        raise SystemExit(f"Terminal edge {int(edges[-1])} != chromosome-size sum {cum}")
    return edges, chromnames, chromsizes, offsets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cooler", required=True, help="Validated raw 250-kb .cool file")
    ap.add_argument("--display-params", required=True, help="Figure 3 parameter/value TSV")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--stem", default="Figure3_code_base_true_coordinates")
    args = ap.parse_args()

    try:
        import cooler
    except ImportError:
        raise SystemExit("Install cooler and h5py before running this script")

    params = read_parameters(Path(args.display_params))
    font_size = float(params.get("font_size_pt", 7.0))
    set_pub_style(font_size)
    if int(params["resolution_bp"]) != 250000:
        raise SystemExit(f"Expected 250000-bp resolution; got {params['resolution_bp']}")
    if params["normalization"].strip().lower() != "none":
        raise SystemExit(f"Expected raw/unbalanced matrix; got {params['normalization']}")
    if params["transform"].strip() != "log10(raw count + 1)":
        raise SystemExit(f"Unexpected transform: {params['transform']}")
    if "true cumulative sequence coordinate" not in params["coordinate_system"].lower():
        raise SystemExit(f"Unexpected coordinate system: {params['coordinate_system']}")

    cool_path = Path(args.cooler)
    expected_sha = params.get("input_cool_sha256", "").strip().lower()
    if expected_sha:
        observed_sha = sha256_file(cool_path)
        if observed_sha.lower() != expected_sha:
            raise SystemExit(
                f"Cooler SHA256 {observed_sha} != declared {expected_sha}"
            )
    c = cooler.Cooler(str(cool_path))
    raw = c.matrix(balance=False)[:]
    z = np.log10(raw.astype(np.float64) + 1.0)
    edges_bp, chromnames, chromsizes, offsets = true_cumulative_edges(c)
    if z.shape != (len(edges_bp)-1, len(edges_bp)-1):
        raise SystemExit(f"Matrix shape {z.shape} does not match {len(edges_bp)-1} bins")

    declared_total = int(params["primary_total_bp"])
    observed_total = int(edges_bp[-1])
    if observed_total != declared_total:
        raise SystemExit(f"Cooler total {observed_total} != declared total {declared_total}")
    total_mb = observed_total / 1e6
    edges_mb = edges_bp / 1e6
    vmin, vmax = float(params["vmin_log10"]), float(params["vmax_log10"])

    width_mm = float(params.get("figure_width_mm", 170.18))
    height_mm = float(params.get("figure_height_mm", 153.67))
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    ax = fig.add_axes([0.105, 0.095, 0.69, 0.79])
    cmap = params.get("colormap", "Reds")
    mesh = ax.pcolormesh(edges_mb, edges_mb, z, cmap=cmap, vmin=vmin, vmax=vmax,
                         shading="flat", rasterized=True)
    ax.set_xlim(0, total_mb)
    ax.set_ylim(total_mb, 0)
    ax.set_aspect("equal", adjustable="box")

    boundaries = [offsets[ch]/1e6 for ch in chromnames] + [total_mb]
    boundary_lw = float(params.get("boundary_linewidth_pt", 0.50))
    for x in boundaries:
        ax.axvline(x, color="black", lw=boundary_lw)
        ax.axhline(x, color="black", lw=boundary_lw)

    centers = [(offsets[ch] + int(chromsizes[ch])/2)/1e6 for ch in chromnames]
    labels = [ch.replace("Chr", "") for ch in chromnames]
    ax.set_xticks(centers, labels)
    ax.set_yticks(centers, labels)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Chromosome")
    ax.tick_params(axis="both", length=0, pad=2)

    coord_ticks = list(np.arange(0, 161, 20, dtype=float))
    if total_mb > coord_ticks[-1] + 0.5:
        coord_ticks.append(total_mb)
    fmt = coord_formatter(total_mb)
    top = ax.secondary_xaxis("top")
    top.set_xticks(coord_ticks)
    top.xaxis.set_major_formatter(fmt)
    top.set_xlabel("Sequence coordinate (Mb)")
    top.tick_params(length=2.2, pad=2)
    right = ax.secondary_yaxis("right")
    right.set_yticks(coord_ticks)
    right.yaxis.set_major_formatter(fmt)
    right.set_ylabel("Sequence coordinate (Mb)")
    right.tick_params(length=2.2, pad=2)

    cbax = fig.add_axes([0.865, 0.11, 0.032, 0.76])
    cb = fig.colorbar(mesh, cax=cbax)
    cb.set_label(r"$\log_{10}(\mathrm{raw\ contact\ count}+1)$", labelpad=6)
    cb.outline.set_linewidth(0.5)
    ticks = np.arange(np.ceil(vmin*2)/2, np.floor(vmax*2)/2 + 0.001, 0.5)
    cb.set_ticks(ticks)
    cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x:.1f}"))

    expected_chroms = [f"Chr{i:02d}" for i in range(1, 15)]
    if chromnames != expected_chroms:
        raise SystemExit(f"Unexpected chromosome order/names: {chromnames}")

    print("Bins:", len(edges_bp)-1)
    print("Total bp:", observed_total)
    print("Display range:", vmin, vmax)
    print("Terminal bin width bp:", int(edges_bp[-1]-edges_bp[-2]))
    print("Figure size mm:", width_mm, height_mm)
    print("Font size pt:", font_size)
    print("Boundary linewidth pt:", boundary_lw)
    if expected_sha:
        print("Cooler SHA256: PASS", expected_sha)
    save_triplet_fixed(fig, args.outdir, args.stem, dpi=int(params.get("png_dpi", 450)))
    plt.close(fig)

if __name__ == "__main__":
    main()

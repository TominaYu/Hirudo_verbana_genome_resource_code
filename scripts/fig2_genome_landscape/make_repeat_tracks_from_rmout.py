#!/usr/bin/env python3
import os
import re
import math

RMOUT = os.environ["RMOUT"]
MAPFILE = os.environ.get("MAPFILE", "data/pri_to_chr.map.tsv")
CHROMS = os.environ.get("CHROMS", "data/chrom.sizes")
OUTDIR = os.environ.get("OUTDIR", "data")
WINDOW = int(os.environ.get("WINDOW", "100000"))

def read_map(path):
    m = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            a, b = line.rstrip("\n").split("\t")
            m[a] = b
    return m

def read_chrom_sizes(path):
    order = []
    sizes = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, size = line.rstrip("\n").split("\t")
            order.append(chrom)
            sizes[chrom] = int(size)
    return order, sizes

def classify_repeat(classfam):
    x = classfam.upper()

    excluded = [
        "SIMPLE_REPEAT", "LOW_COMPLEXITY", "SATELLITE",
        "MICROSATELLITE", "TANDEM", "TRF",
        "RRNA", "TRNA", "SNRNA", "SCRNA", "SRPRNA"
    ]
    if any(e in x for e in excluded):
        return None

    out = set()

    if "LINE" in x:
        out.add("LINE")
    if "LTR" in x:
        out.add("LTR")

    te_keywords = [
        "LINE", "LTR", "SINE", "DNA", "RC", "HELITRON",
        "ROLLING", "PENELOPE", "RETROPOSON", "GYPSY",
        "COPIA", "BEL", "DIRS", "ERV"
    ]
    if any(k in x for k in te_keywords):
        out.add("TE")

    return out if out else None

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]

def write_track(path, chroms, sizes, values):
    with open(path, "w") as out:
        for chrom in chroms:
            for i, val in enumerate(values[chrom]):
                s = i * WINDOW
                e = min((i + 1) * WINDOW, sizes[chrom])
                out.write(f"{chrom}\t{s}\t{e}\t{val:.6f}\n")

def max_track(track):
    return max(v for chrom in track for v in track[chrom])

def main():
    name_map = read_map(MAPFILE)
    chroms, sizes = read_chrom_sizes(CHROMS)
    chromset = set(chroms)

    rep = {track: {c: [] for c in chroms} for track in ["LINE", "LTR", "TE"]}

    with open(RMOUT) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("SW") or line.startswith("score") or line.startswith("There were"):
                continue
            if line.startswith("#"):
                continue

            parts = re.split(r"\s+", line.strip())
            if len(parts) < 11:
                continue

            query = parts[4]
            if query not in name_map:
                continue

            chrom = name_map[query]
            if chrom not in chromset:
                continue

            try:
                start = int(parts[5]) - 1
                end = int(parts[6])
            except ValueError:
                continue

            classfam = parts[10]
            classes = classify_repeat(classfam)
            if not classes:
                continue

            for cls in classes:
                rep[cls][chrom].append((start, end))

    line_frac = {}
    ltr_frac = {}
    te_frac = {}

    for chrom in chroms:
        size = sizes[chrom]
        nwin = math.ceil(size / WINDOW)

        line_bp = [0] * nwin
        ltr_bp = [0] * nwin
        te_bp = [0] * nwin

        merged_line = merge_intervals(rep["LINE"][chrom])
        merged_ltr  = merge_intervals(rep["LTR"][chrom])
        merged_te   = merge_intervals(rep["TE"][chrom])

        def add_coverage(intervals, arr):
            for s, e in intervals:
                w0 = s // WINDOW
                w1 = (e - 1) // WINDOW
                for w in range(w0, w1 + 1):
                    ws = w * WINDOW
                    we = min((w + 1) * WINDOW, size)
                    ov = max(0, min(e, we) - max(s, ws))
                    arr[w] += ov

        add_coverage(merged_line, line_bp)
        add_coverage(merged_ltr,  ltr_bp)
        add_coverage(merged_te,   te_bp)

        line_frac[chrom] = []
        ltr_frac[chrom] = []
        te_frac[chrom] = []

        for w in range(nwin):
            ws = w * WINDOW
            we = min((w + 1) * WINDOW, size)
            width = we - ws
            line_frac[chrom].append(line_bp[w] / width if width > 0 else 0.0)
            ltr_frac[chrom].append(ltr_bp[w] / width if width > 0 else 0.0)
            te_frac[chrom].append(te_bp[w] / width if width > 0 else 0.0)

    write_track(f"{OUTDIR}/line_fraction.txt", chroms, sizes, line_frac)
    write_track(f"{OUTDIR}/ltr_fraction.txt", chroms, sizes, ltr_frac)
    write_track(f"{OUTDIR}/te_fraction.txt", chroms, sizes, te_frac)

    with open(f"{OUTDIR}/repeat_summary.txt", "w") as out:
        out.write(f"LINEfrac_max\t{max_track(line_frac):.6f}\n")
        out.write(f"LTRfrac_max\t{max_track(ltr_frac):.6f}\n")
        out.write(f"TEfrac_max\t{max_track(te_frac):.6f}\n")

if __name__ == "__main__":
    main()

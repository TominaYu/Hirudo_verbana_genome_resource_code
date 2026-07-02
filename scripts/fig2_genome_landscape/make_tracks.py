#!/usr/bin/env python3
import os
import re
import math

FA = os.environ["FA"]
GTF = os.environ["GTF"]
RMGFF = os.environ["RMGFF"]
CHROMS = os.environ.get("CHROMS", "data/chrom.sizes")
OUTDIR = os.environ.get("OUTDIR", "data")
WINDOW = int(os.environ.get("WINDOW", "100000"))

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

def read_fasta_subset(path, chrom_order):
    keep = set(chrom_order)
    seqs = {}
    name = None
    buf = []
    use = False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if use and name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].strip().split()[0]
                use = name in keep
                buf = []
            else:
                if use:
                    buf.append(line.strip())
        if use and name is not None:
            seqs[name] = "".join(buf).upper()
    missing = [c for c in chrom_order if c not in seqs]
    if missing:
        raise RuntimeError(f"FASTA に見つからない染色体: {missing}")
    return seqs

def parse_gene_id(attr):
    pats = [
        r'gene_id "([^"]+)"',
        r'gene_id=([^;]+)',
        r'ID=([^;]+)',
        r'Parent=([^;]+)'
    ]
    for p in pats:
        m = re.search(p, attr)
        if m:
            return m.group(1)
    return None

def parse_gtf_gene_spans(path, chroms):
    chromset = set(chroms)
    spans = {c: {} for c in chroms}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, _, _, attr = fields
            if chrom not in chromset:
                continue
            if feature not in {"gene", "transcript", "mRNA", "exon", "CDS"}:
                continue
            gid = parse_gene_id(attr)
            if gid is None:
                continue
            s = int(start) - 1
            e = int(end)
            if gid not in spans[chrom]:
                spans[chrom][gid] = [s, e]
            else:
                if s < spans[chrom][gid][0]:
                    spans[chrom][gid][0] = s
                if e > spans[chrom][gid][1]:
                    spans[chrom][gid][1] = e
    return spans

def classify_repeat(attr, source="", feature=""):
    text = f"{source};{feature};{attr}".upper()
    text = text.replace(" ", "_")

    excluded = [
        "SIMPLE_REPEAT", "LOW_COMPLEXITY", "SATELLITE",
        "MICROSATELLITE", "TANDEM_REPEAT", "TRF"
    ]
    te_patterns = [
        r'(^|[;:/,_-])LINE([;:/,_-]|$)',
        r'(^|[;:/,_-])LTR([;:/,_-]|$)',
        r'(^|[;:/,_-])SINE([;:/,_-]|$)',
        r'(^|[;:/,_-])DNA([;:/,_-]|$)',
        r'HELITRON',
        r'ROLLING_CIRCLE',
        r'PENELOPE',
        r'RETROPOSON',
        r'GYPSY',
        r'COPIA',
        r'BEL',
        r'DIRS',
        r'ERV'
    ]

    is_line = re.search(r'(^|[;:/,_-])LINE([;:/,_-]|$)', text) is not None
    is_ltr  = re.search(r'(^|[;:/,_-])LTR([;:/,_-]|$)', text) is not None
    is_te   = any(re.search(p, text) for p in te_patterns)

    if any(x in text for x in excluded) and not (is_line or is_ltr or is_te):
        return None

    out = set()
    if is_line:
        out.add("LINE")
    if is_ltr:
        out.add("LTR")
    if is_te or is_line or is_ltr:
        out.add("TE")
    return out if out else None

def read_repeat_intervals(path, chroms):
    chromset = set(chroms)
    rep = {track: {c: [] for c in chroms} for track in ["LINE", "LTR", "TE"]}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, _, _, _, attr = fields
            if chrom not in chromset:
                continue
            classes = classify_repeat(attr, source=source, feature=feature)
            if not classes:
                continue
            s = int(start) - 1
            e = int(end)
            for cls in classes:
                rep[cls][chrom].append((s, e))
    return rep

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

def make_windows(size, window):
    n = math.ceil(size / window)
    return [(i * window, min((i + 1) * window, size)) for i in range(n)]

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
    chroms, sizes = read_chrom_sizes(CHROMS)
    seqs = read_fasta_subset(FA, chroms)

    gc_skew = {}
    gc_ratio = {}
    for chrom in chroms:
        size = sizes[chrom]
        seq = seqs[chrom][:size]
        windows = make_windows(size, WINDOW)
        gc_skew[chrom] = []
        gc_ratio[chrom] = []
        for s, e in windows:
            sub = seq[s:e]
            a = sub.count("A")
            t = sub.count("T")
            g = sub.count("G")
            c = sub.count("C")
            denom_gc = g + c
            denom_all = a + t + g + c
            skew = (g - c) / denom_gc if denom_gc > 0 else 0.0
            ratio = denom_gc / denom_all if denom_all > 0 else 0.0
            gc_skew[chrom].append(skew)
            gc_ratio[chrom].append(ratio)

    gene_spans = parse_gtf_gene_spans(GTF, chroms)
    gene_density = {}
    total_genes = 0
    for chrom in chroms:
        nwin = math.ceil(sizes[chrom] / WINDOW)
        gene_density[chrom] = [0] * nwin
        for gid, (s, e) in gene_spans[chrom].items():
            total_genes += 1
            w0 = s // WINDOW
            w1 = (e - 1) // WINDOW
            for w in range(w0, w1 + 1):
                gene_density[chrom][w] += 1

    rep = read_repeat_intervals(RMGFF, chroms)
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

    write_track(f"{OUTDIR}/gc_skew.txt", chroms, sizes, gc_skew)
    write_track(f"{OUTDIR}/gc_ratio.txt", chroms, sizes, gc_ratio)
    write_track(f"{OUTDIR}/gene_density.txt", chroms, sizes, gene_density)
    write_track(f"{OUTDIR}/line_fraction.txt", chroms, sizes, line_frac)
    write_track(f"{OUTDIR}/ltr_fraction.txt", chroms, sizes, ltr_frac)
    write_track(f"{OUTDIR}/te_fraction.txt", chroms, sizes, te_frac)

    colors = [
        "chr01","chr02","chr03","chr04","chr05","chr06","chr07",
        "chr08","chr09","chr10","chr11","chr12","chr13","chr14"
    ]
    with open(f"{OUTDIR}/karyotype.hirudo.txt", "w") as out:
        for chrom, color in zip(chroms, colors):
            out.write(f"chr - {chrom} {chrom} 0 {sizes[chrom]} {color}\n")

    with open(f"{OUTDIR}/summary.txt", "w") as out:
        out.write(f"Genes_total\t{total_genes}\n")
        out.write(f"GCratio_max\t{max_track(gc_ratio):.6f}\n")
        out.write(f"GCskew_max\t{max_track(gc_skew):.6f}\n")
        out.write(f"GeneDensity_max\t{max_track(gene_density):.6f}\n")
        out.write(f"LINEfrac_max\t{max_track(line_frac):.6f}\n")
        out.write(f"LTRfrac_max\t{max_track(ltr_frac):.6f}\n")
        out.write(f"TEfrac_max\t{max_track(te_frac):.6f}\n")

if __name__ == "__main__":
    main()

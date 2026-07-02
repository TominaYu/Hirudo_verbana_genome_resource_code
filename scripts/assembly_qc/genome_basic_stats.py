#!/usr/bin/env python3
import argparse
from collections import OrderedDict

def read_fasta(path):
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, "".join(seq)

def n50_l50(lengths):
    if not lengths:
        return 0, 0
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    half = total / 2
    acc = 0
    for i, x in enumerate(lengths, start=1):
        acc += x
        if acc >= half:
            return x, i
    return 0, 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-lengths-tsv", required=True)
    args = ap.parse_args()

    lengths = []
    total_len = 0
    total_gc = 0
    total_acgt = 0
    total_n = 0
    total_other = 0
    total_softmasked = 0

    per_record = []

    for name, seq in read_fasta(args.fasta):
        L = len(seq)
        lengths.append(L)
        total_len += L

        seq_u = seq.upper()
        a = seq_u.count("A")
        c = seq_u.count("C")
        g = seq_u.count("G")
        t = seq_u.count("T")
        n = seq_u.count("N")
        acgt = a + c + g + t
        other = L - acgt - n
        gc = g + c
        softmasked = sum(1 for ch in seq if ch.islower())

        total_gc += gc
        total_acgt += acgt
        total_n += n
        total_other += other
        total_softmasked += softmasked

        per_record.append((name, L, gc, acgt, n, other, softmasked))

    record_count = len(lengths)
    n50, l50 = n50_l50(lengths)
    longest = max(lengths) if lengths else 0
    shortest = min(lengths) if lengths else 0

    gc_pct_acgt = (100.0 * total_gc / total_acgt) if total_acgt else 0.0
    gc_pct_total = (100.0 * total_gc / total_len) if total_len else 0.0
    n_pct_total = (100.0 * total_n / total_len) if total_len else 0.0
    other_pct_total = (100.0 * total_other / total_len) if total_len else 0.0
    softmasked_pct_total = (100.0 * total_softmasked / total_len) if total_len else 0.0

    with open(args.out_tsv, "w") as f:
        f.write("label\trecord_count\ttotal_bp\tgc_pct_acgt\tgc_pct_total\tn_pct_total\tother_pct_total\tsoftmasked_pct_total\tn50\tl50\tlongest\tshortest\n")
        f.write(
            f"{args.label}\t{record_count}\t{total_len}\t"
            f"{gc_pct_acgt:.6f}\t{gc_pct_total:.6f}\t{n_pct_total:.6f}\t{other_pct_total:.6f}\t{softmasked_pct_total:.6f}\t"
            f"{n50}\t{l50}\t{longest}\t{shortest}\n"
        )

    with open(args.out_lengths_tsv, "w") as f:
        f.write("record_id\tlength_bp\tgc_count\tacgt_count\tn_count\tother_count\tsoftmasked_count\n")
        for row in per_record:
            f.write("\t".join(map(str, row)) + "\n")

if __name__ == "__main__":
    main()

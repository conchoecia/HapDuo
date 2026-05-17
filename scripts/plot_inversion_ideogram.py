#!/usr/bin/env python3
"""
Inversion ideogram: 14 hapA chromosomes as horizontal bars, with all called
inversions overlaid as colored rectangles at their hapA t-coordinates.

Color encodes inversion size (>=1 Mb, 100 kb-1 Mb, <100 kb).
A small inset histogram shows the per-Mb inversion density per chromosome.
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

CHR_RE = re.compile(r"_chr(\d+)$")

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'


def chr_num(name):
    m = CHR_RE.search(name)
    return int(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--breakpoints", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--hap-prefix", default="eupHapAv0.3_chr")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-anchors", type=int, default=1,
                    help="ignore inversions supported by fewer than this many "
                         "chained antiparallel anchors. 3 = excludes single-"
                         "block calls (Hi-C-visible floor).")
    ap.add_argument("--min-len", type=int, default=0,
                    help="ignore inversions below this length in bp")
    args = ap.parse_args()

    # chromosome lengths from fai
    chr_len = {}
    for ln in open(args.fai):
        name, length = ln.split("\t")[0], int(ln.split("\t")[1])
        if name.startswith(args.hap_prefix):
            chr_len[name] = length
    chrs = sorted(chr_len.keys(), key=lambda c: chr_num(c))

    # load inversions (left+right rows per rank)
    pairs = defaultdict(dict)
    with open(args.breakpoints) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pairs[int(row["rank"])][row["end"]] = row
    inversions = []
    for rank, sides in pairs.items():
        L, R = sides["left"], sides["right"]
        chr_a = L["hapA_chr"]
        t_lo = min(int(L["hapA_pos"]), int(R["hapA_pos"]))
        t_hi = max(int(L["hapA_pos"]), int(R["hapA_pos"]))
        inv_len = int(L["inv_len"])
        n_anch = int(L["n_anchors"])
        if n_anch < args.min_anchors or inv_len < args.min_len:
            continue
        inversions.append(dict(rank=rank, chr=chr_a, t_lo=t_lo, t_hi=t_hi,
                               inv_len=inv_len, n_anchors=n_anch))

    # size bins for color
    def color_for(inv_len):
        if inv_len >= 1_000_000:    return "#b30000"  # large red
        if inv_len >=   100_000:    return "#fb6a4a"  # mid orange
        return "#fcae91"                              # small pale

    # counts per chr per size class
    per_chr = defaultdict(lambda: dict(big=0, mid=0, small=0))
    for inv in inversions:
        if inv["inv_len"] >= 1_000_000:
            per_chr[inv["chr"]]["big"] += 1
        elif inv["inv_len"] >= 100_000:
            per_chr[inv["chr"]]["mid"] += 1
        else:
            per_chr[inv["chr"]]["small"] += 1

    # ---------- figure ----------
    fig = plt.figure(figsize=(8.5, 6.5))
    ax = fig.add_axes([0.10, 0.10, 0.78, 0.85])
    inset = fig.add_axes([0.90, 0.10, 0.08, 0.85])

    max_len = max(chr_len[c] for c in chrs)
    chr_h = 0.55
    bp_to_mb = 1e-6

    for i, c in enumerate(chrs):
        y = len(chrs) - 1 - i  # chr1 at top
        L = chr_len[c]
        # chromosome backbone
        ax.add_patch(mpatches.FancyBboxPatch(
            (0, y - chr_h / 2), L * bp_to_mb, chr_h,
            boxstyle="round,pad=0,rounding_size=0.2",
            facecolor="#eeeeee", edgecolor="#777777", linewidth=0.7))
        # inversions on this chr
        my_invs = [iv for iv in inversions if iv["chr"] == c]
        my_invs.sort(key=lambda iv: -iv["inv_len"])  # big behind small? draw big first
        for iv in my_invs:
            x_lo = iv["t_lo"] * bp_to_mb
            x_hi = iv["t_hi"] * bp_to_mb
            ax.add_patch(mpatches.Rectangle(
                (x_lo, y - chr_h / 2), max(0.05, x_hi - x_lo), chr_h,
                facecolor=color_for(iv["inv_len"]),
                edgecolor="#5a0000", linewidth=0.3, alpha=0.85))
        # chr label on the left
        ax.text(-0.5, y, c.replace(args.hap_prefix, "chr"),
                ha="right", va="center", fontsize=9, color="#333333")
        # length label on the right
        ax.text(L * bp_to_mb + 0.3, y, f"{L/1e6:.1f} Mb",
                ha="left", va="center", fontsize=7, color="#888888")

    ax.set_xlim(-1, max_len * bp_to_mb + 4)
    ax.set_ylim(-1, len(chrs))
    ax.set_xlabel("haplotype-A position (Mb)")
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#555555")

    # legend (top-left over empty space)
    legend = [
        mpatches.Patch(color="#b30000", label="≥ 1 Mb"),
        mpatches.Patch(color="#fb6a4a", label="100 kb – 1 Mb"),
        mpatches.Patch(color="#fcae91", label="< 100 kb"),
    ]
    ax.legend(handles=legend, loc="upper right", fontsize=8,
              frameon=True, framealpha=0.9, edgecolor="#bbbbbb",
              title="inversion size", title_fontsize=8)

    # inset: counts per chr
    counts_big = [per_chr[c]["big"] for c in chrs]
    counts_mid = [per_chr[c]["mid"] for c in chrs]
    counts_sm  = [per_chr[c]["small"] for c in chrs]
    ys = np.arange(len(chrs))[::-1]
    inset.barh(ys, counts_big, color="#b30000", height=0.7,
               label="≥1 Mb")
    inset.barh(ys, counts_mid, left=counts_big, color="#fb6a4a",
               height=0.7, label="100k–1M")
    inset.barh(ys, counts_sm,
               left=np.array(counts_big) + np.array(counts_mid),
               color="#fcae91", height=0.7, label="<100k")
    inset.set_yticks([])
    inset.set_ylim(-1, len(chrs))
    inset.set_xlabel("count", fontsize=7)
    inset.tick_params(axis="x", labelsize=6)
    for s in ("top", "right"):
        inset.spines[s].set_visible(False)

    # summary in caption space
    n_total = len(inversions)
    n_big = sum(1 for iv in inversions if iv["inv_len"] >= 1_000_000)
    n_mid = sum(1 for iv in inversions
                if 100_000 <= iv["inv_len"] < 1_000_000)
    n_sm  = n_total - n_big - n_mid
    total_inv_bp = sum(iv["inv_len"] for iv in inversions)
    total_chr_bp = sum(chr_len[c] for c in chrs)
    filter_txt = ""
    if args.min_anchors > 1:
        filter_txt += f", ≥{args.min_anchors} antiparallel anchors"
    if args.min_len > 0:
        filter_txt += f", ≥{args.min_len/1000:g} kb"
    fig.suptitle(
        f"High-confidence inter-haplotype inversions (n={n_total}{filter_txt}): "
        f"{n_big} ≥1 Mb, {n_mid} 100 kb–1 Mb, {n_sm} <100 kb  |  "
        f"total inverted span {total_inv_bp/1e6:.1f} Mb "
        f"({total_inv_bp/total_chr_bp*100:.1f}% of hapA chr length)",
        fontsize=9, y=0.98)

    out_pdf = args.out_prefix + ".pdf"
    out_png = args.out_prefix + ".png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    print(f"# wrote {out_pdf}", file=sys.stderr)
    print(f"# wrote {out_png}", file=sys.stderr)


if __name__ == "__main__":
    main()

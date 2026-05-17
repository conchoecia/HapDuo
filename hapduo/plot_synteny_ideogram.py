#!/usr/bin/env python3
"""
Per-chromosome two-haplotype synteny ideogram.

Same row-per-chromosome layout as the inversion-call ideogram, but each row
shows both the haplotype-A chromosome bar (top) and the haplotype-B
chromosome bar (bottom), drawn to their true Mb length on a shared x-axis,
with every PAF alignment block overlaid as a coloured polygon ribbon
connecting the two bars.

  - Ribbons are coloured by strand: blue for `+` (colinear), red for `-`
    (antiparallel / inverted).
  - `--min-block` filters out tiny noise alignments so the figure is
    legible at chromosome scale.
  - PAF records on mis-matched chromosome pairs (e.g.\ hapA chr2 → hapB
    chr11) are not drawn in this layout because each row is a chr-pair;
    the count of such cross-chromosome blocks is reported in the title.

Usage:
  plot_synteny_ideogram.py
    --paf  hapB_to_hapA.paf
    --fai  eup_hapAhapBMetaExtraMito_v0.3.fasta.fai
    --out-prefix out/synteny_ideogram
    --min-block 50000
    --chain
"""
from __future__ import annotations
import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PolyCollection
import numpy as np

from .paf_chain import parse_paf, chain_paf

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"

CHR_RE   = re.compile(r"_chr(\d+)$")
FWD_COLOR = "#3a6fa8"
REV_COLOR = "#c2484d"
CHR_FILL  = "#eeeeee"
CHR_EDGE  = "#666666"


def chr_num(name):
    m = CHR_RE.search(name)
    return int(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--hapA-prefix", default="eupHapAv0.3_chr")
    ap.add_argument("--hapB-prefix", default="eupHapBv0.3_chr")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-block", type=int, default=50_000,
                    help="ignore PAF records shorter than this on either axis")
    ap.add_argument("--chain", action="store_true",
                    help="run the same greedy chainer the inversion detector "
                         "uses, so a long syntenic block is drawn as one ribbon")
    ap.add_argument("--alpha", type=float, default=0.55)
    ap.add_argument("--figsize", nargs=2, type=float, default=[10, 8])
    args = ap.parse_args()

    # --- chromosome lengths from fai ---
    chr_lens_a, chr_lens_b = {}, {}
    for ln in open(args.fai):
        name, length = ln.split("\t")[0], int(ln.split("\t")[1])
        if name.startswith(args.hapA_prefix):
            chr_lens_a[name] = length
        elif name.startswith(args.hapB_prefix):
            chr_lens_b[name] = length
    # match by chr number
    chr_nums = sorted({chr_num(c) for c in chr_lens_a} | {chr_num(c) for c in chr_lens_b})
    print(f"# {len(chr_nums)} chromosome pairs", file=sys.stderr)
    pairs = []
    for n in chr_nums:
        a = [c for c in chr_lens_a if chr_num(c) == n]
        b = [c for c in chr_lens_b if chr_num(c) == n]
        if a and b:
            pairs.append((n, a[0], b[0]))

    # max chromosome length across both haplotypes — sets the x-axis range
    max_len = max(max(chr_lens_a.values()), max(chr_lens_b.values()))
    bp_to_mb = 1e-6

    # --- PAF ---
    blocks = parse_paf(args.paf)
    print(f"# {len(blocks)} PAF records loaded", file=sys.stderr)
    if args.chain:
        blocks = chain_paf(blocks, max_gap=50_000, rel_skew=0.3,
                           min_block_len=1000)
        print(f"# {len(blocks)} chained syntenic blocks", file=sys.stderr)

    # bucket blocks by chr-pair
    by_pair = defaultdict(list)
    n_kept = n_trans = 0
    for b in blocks:
        t_len = b.t_end - b.t_start
        q_len = b.q_end - b.q_start
        if min(t_len, q_len) < args.min_block:
            continue
        if b.tname not in chr_lens_a or b.qname not in chr_lens_b:
            continue
        if chr_num(b.tname) != chr_num(b.qname):
            n_trans += 1
            continue
        by_pair[chr_num(b.tname)].append(b)
        n_kept += 1
    print(f"# {n_kept} blocks kept (>= {args.min_block:,} bp on shorter axis), "
          f"{n_trans} cross-chromosome blocks skipped", file=sys.stderr)

    # ---------------- layout ----------------
    n_rows = len(pairs)
    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.94, bottom=0.05)

    row_h    = 1.0         # row spacing
    chr_h    = 0.22        # height of each chromosome bar within a row
    a_off    = 0.30        # vertical offset of hapA bar centre above the row centre
    b_off    = -0.30       # vertical offset of hapB bar centre below the row centre

    n_fwd_drawn = n_rev_drawn = 0

    for i, (n, a_name, b_name) in enumerate(pairs):
        y_row = (n_rows - 1 - i) * row_h          # chr1 at top
        ya = y_row + a_off
        yb = y_row + b_off

        # hapA chromosome bar
        wa = chr_lens_a[a_name] * bp_to_mb
        ax.add_patch(mpatches.Rectangle(
            (0, ya - chr_h / 2), wa, chr_h,
            facecolor=CHR_FILL, edgecolor=CHR_EDGE, lw=0.5))

        # hapB chromosome bar
        wb = chr_lens_b[b_name] * bp_to_mb
        ax.add_patch(mpatches.Rectangle(
            (0, yb - chr_h / 2), wb, chr_h,
            facecolor=CHR_FILL, edgecolor=CHR_EDGE, lw=0.5))

        # left-side row label
        ax.text(-1.5, y_row, f"chr{n}", ha="right", va="center",
                fontsize=9, color="#333333")

        # tiny haplotype labels on the left of each bar
        ax.text(-0.3, ya, "A", ha="right", va="center",
                fontsize=7, color="#1f4f8a")
        ax.text(-0.3, yb, "B", ha="right", va="center",
                fontsize=7, color="#226a1f")

        # ribbons
        fwd_polys, rev_polys = [], []
        for b in by_pair.get(n, []):
            xa1 = b.t_start * bp_to_mb
            xa2 = b.t_end   * bp_to_mb
            if b.strand == "+":
                xb1 = b.q_start * bp_to_mb
                xb2 = b.q_end   * bp_to_mb
            else:
                xb1 = b.q_end   * bp_to_mb
                xb2 = b.q_start * bp_to_mb
            poly = [(xa1, ya - chr_h / 2),
                    (xa2, ya - chr_h / 2),
                    (xb2, yb + chr_h / 2),
                    (xb1, yb + chr_h / 2)]
            if b.strand == "+":
                fwd_polys.append(poly)
                n_fwd_drawn += 1
            else:
                rev_polys.append(poly)
                n_rev_drawn += 1
        if fwd_polys:
            ax.add_collection(PolyCollection(fwd_polys, facecolors=FWD_COLOR,
                                             edgecolors="none", alpha=args.alpha,
                                             antialiased=False, rasterized=True))
        if rev_polys:
            ax.add_collection(PolyCollection(rev_polys, facecolors=REV_COLOR,
                                             edgecolors="none", alpha=args.alpha,
                                             antialiased=False, rasterized=True))

    # x-axis: shared Mb scale
    ax.set_xlim(-3, max_len * bp_to_mb + 0.5)
    ax.set_ylim(-1.2, n_rows * row_h - 0.4)
    ax.set_xlabel("chromosome position (Mb)", fontsize=9)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#555555")
    ax.tick_params(axis="x", labelsize=8)

    # legend (top-right of the figure)
    fig.legend(handles=[
        mpatches.Patch(color=FWD_COLOR, label="alignment (+) colinear"),
        mpatches.Patch(color=REV_COLOR, label="alignment (-) inverted"),
    ], loc="upper right", bbox_to_anchor=(0.99, 0.985), fontsize=8,
       frameon=True, framealpha=0.9, edgecolor="#bbbbbb")

    # title summary
    fig.suptitle(
        f"Two-haplotype synteny per chromosome  "
        f"({n_kept:,} blocks $\\geq${args.min_block/1000:g} kb: "
        f"{n_fwd_drawn:,} colinear, {n_rev_drawn:,} inverted)",
        fontsize=10, y=0.985, x=0.08, ha="left")

    out_pdf = Path(args.out_prefix + ".pdf")
    out_png = Path(args.out_prefix + ".png")
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=220)
    print(f"# wrote {out_pdf}", file=sys.stderr)
    print(f"# wrote {out_png}", file=sys.stderr)


if __name__ == "__main__":
    main()

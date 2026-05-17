#!/usr/bin/env python3
"""
Style 4: full-page inversion plot with Hi-C support.

Top → bottom layout:
  hapA Hi-C pyramid heatmap  (white→red, log-counts)
  hapA CCS read pile
  hapA chromosome track   (with vertical breakpoint markers)
  ribbons of chained PAF alignments (blue = colinear, red = inverted)
  hapB chromosome track   (with vertical breakpoint markers)
  hapB CCS read pile
  hapB Hi-C pyramid heatmap  (mirrored upside-down, white→red)

The Hi-C file from a competitive-mapping run treats the entire concatenated
assembly as one pseudo-contig named "assembly" or similar — pass the
concatenation offsets of the hapA and hapB scaffolds with --hic-offset-a /
--hic-offset-b.

Usage example:
  plot_style4_hic.py
    --paf hapB_to_hapA.paf
    --ccs CCS.bam
    --fai REF.fasta.fai
    --hic euplov03_allLibs_q_5.hic
    --hic-chrom assembly
    --hic-offset-a 52815143    # hapA chr offset within concat ref
    --hic-offset-b 243630643   # hapB chr offset within concat ref
    --hic-resolution 10000
    --hapA-scaf eupHapAv0.3_chr4 --hapA-left 3905869 --hapA-right 4229847
    --hapB-scaf eupHapBv0.3_chr4 --hapB-left 6500377 --hapB-right 6207152
    --out-prefix out/style4/rank3_chr4
"""
from __future__ import annotations
import argparse
import os
import sys

from .paf_chain import parse_paf, chain_paf

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.collections import PolyCollection
import numpy as np
import pysam
import hicstraw

matplotlib.rcParams['pdf.fonttype'] = 42   # TrueType -> editable text in Illustrator
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family']  = 'sans-serif'


FWD_COLOR = "#5079a8"
REV_COLOR = "#c2484d"
READ_COLOR_BG = "#e2e2e2"
READ_COLOR_CROSS = "#222222"
BREAK_COLOR = "#b53630"
TRACK_COLOR = "#444444"
WHITE_RED = LinearSegmentedColormap.from_list("wr", ["white", "#b30000"])


def load_chained(paf, tscaf, qscaf, max_gap, rel_skew, min_input_block,
                 min_output_block):
    all_blocks = parse_paf(paf)
    rel = [b for b in all_blocks if b.tname == tscaf and b.qname == qscaf]
    rel = chain_paf(rel, max_gap=max_gap, rel_skew=rel_skew,
                    min_block_len=min_input_block)
    return [b for b in rel if b.blocklen >= min_output_block]


def draw_chr_cartoon(ax, paf_path, hapA_scaf, hapB_scaf, hapA_len, hapB_len,
                     hapA_left, hapA_right, hapB_left, hapB_right,
                     min_block=5000, chain=True):
    """Render a tiny two-haplotype synteny strip for one chromosome with the
    current inversion outlined. Drawn into the supplied inset Axes."""
    blks = parse_paf(paf_path)
    blks = [b for b in blks if b.tname == hapA_scaf and b.qname == hapB_scaf]
    if chain:
        blks = chain_paf(blks, max_gap=50_000, rel_skew=0.3,
                         min_block_len=1000)
    blks = [b for b in blks
            if min(b.t_end - b.t_start, b.q_end - b.q_start) >= min_block]

    # Layout: hapA bar at y=0.75, hapB bar at y=0.25; bars sized by chr length
    # against a shared x-scale [0, max_len]
    max_len = max(hapA_len, hapB_len)
    chr_h = 0.18
    A_Y, B_Y = 0.78, 0.22

    def xa(p): return p / max_len
    def xb(p): return p / max_len

    # chromosome backbones
    ax.add_patch(mpatches.Rectangle((0, A_Y - chr_h / 2), hapA_len / max_len,
                                    chr_h, facecolor="#eeeeee",
                                    edgecolor="#666666", lw=0.4))
    ax.add_patch(mpatches.Rectangle((0, B_Y - chr_h / 2), hapB_len / max_len,
                                    chr_h, facecolor="#eeeeee",
                                    edgecolor="#666666", lw=0.4))

    # ribbons
    fwd_polys, rev_polys = [], []
    for b in blks:
        ax1 = xa(b.t_start); ax2 = xa(b.t_end)
        if b.strand == "+":
            bx1 = xb(b.q_start); bx2 = xb(b.q_end)
            fwd_polys.append([(ax1, A_Y - chr_h / 2), (ax2, A_Y - chr_h / 2),
                              (bx2, B_Y + chr_h / 2), (bx1, B_Y + chr_h / 2)])
        else:
            bx1 = xb(b.q_end); bx2 = xb(b.q_start)
            rev_polys.append([(ax1, A_Y - chr_h / 2), (ax2, A_Y - chr_h / 2),
                              (bx2, B_Y + chr_h / 2), (bx1, B_Y + chr_h / 2)])
    if fwd_polys:
        ax.add_collection(PolyCollection(fwd_polys, facecolors=FWD_COLOR,
                                         edgecolors="none", alpha=0.55,
                                         antialiased=False, rasterized=True))
    if rev_polys:
        ax.add_collection(PolyCollection(rev_polys, facecolors=REV_COLOR,
                                         edgecolors="none", alpha=0.55,
                                         antialiased=False, rasterized=True))

    # CALLOUT: outline the current inversion on both bars with a thick black
    # rectangle so the reader can see where in the chromosome they are.
    a_lo, a_hi = sorted([hapA_left, hapA_right])
    b_lo, b_hi = sorted([hapB_left, hapB_right])
    ax.add_patch(mpatches.Rectangle((a_lo / max_len, A_Y - chr_h / 2 - 0.04),
                                    (a_hi - a_lo) / max_len, chr_h + 0.08,
                                    facecolor="none", edgecolor="#000000",
                                    lw=1.2))
    ax.add_patch(mpatches.Rectangle((b_lo / max_len, B_Y - chr_h / 2 - 0.04),
                                    (b_hi - b_lo) / max_len, chr_h + 0.08,
                                    facecolor="none", edgecolor="#000000",
                                    lw=1.2))

    # tiny haplotype labels
    ax.text(-0.02, A_Y, "A", ha="right", va="center", fontsize=6,
            color="#1f4f8a", fontweight="bold")
    ax.text(-0.02, B_Y, "B", ha="right", va="center", fontsize=6,
            color="#226a1f", fontweight="bold")
    ax.text(0.5, 1.01, f"{hapA_scaf.replace('eupHapAv0.3_', '')} (inversion location)",
            ha="center", va="bottom", fontsize=6.5,
            color="#444444", transform=ax.transAxes)

    ax.set_xlim(-0.03, 1.01)
    ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def stack_reads(reads):
    reads = sorted(reads, key=lambda r: r["start"])
    rows = []
    placed = []
    for r in reads:
        for ri, end in enumerate(rows):
            if r["start"] > end:
                rows[ri] = r["end"]
                placed.append((r, ri))
                break
        else:
            rows.append(r["end"])
            placed.append((r, len(rows) - 1))
    return placed, len(rows)


def fetch_coverage(bam, scaf, axis_min, axis_max, n_bins):
    """Per-bin mean CCS read depth across [axis_min, axis_max], counting only
    primary mappings. Returns (bin_centers, depth_per_bin)."""
    bin_width = max(1, (axis_max - axis_min) // n_bins)
    n = max(1, (axis_max - axis_min) // bin_width)
    cov = np.zeros(n, dtype=np.int64)
    for read in bam.fetch(scaf, axis_min, axis_max):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        s = max(axis_min, read.reference_start) - axis_min
        e = min(axis_max, read.reference_end) - axis_min
        if e <= s:
            continue
        i0 = s // bin_width
        i1 = min(n, (e - 1) // bin_width + 1)
        cov[i0:i1] += 1
    centers = axis_min + (np.arange(n) + 0.5) * bin_width
    return centers, cov.astype(float), bin_width


def fetch_reads(bam, scaf, axis_min, axis_max, break_positions,
                break_window=1000):
    """A read is flagged as 'crosses a breakpoint' if its mapped span
    overlaps the (break_window) bp interval centered on any breakpoint."""
    half = break_window // 2
    crossers = set()
    for bp in break_positions:
        for read in bam.fetch(scaf, max(0, bp - half), bp + half):
            crossers.add(read.query_name)
    out = []
    for read in bam.fetch(scaf, axis_min, axis_max):
        out.append({"name": read.query_name,
                    "start": read.reference_start,
                    "end":   read.reference_end,
                    "cross": read.query_name in crossers})
    return out


def fetch_hic_pyramid(hicfile, chrom, offset, start, end, resolution):
    """Return (matrix, bin_starts) for cis contacts in [start, end] at the
    given resolution. Matrix is square len(bin_starts) x len(bin_starts).

    `offset` is the displacement of the local scaffold within the .hic's
    pseudo-contig coordinate system.
    """
    a_start = offset + start
    a_end   = offset + end
    mzd = hicfile.getMatrixZoomData(chrom, chrom, "observed", "NONE",
                                    "BP", resolution)
    arr = mzd.getRecordsAsMatrix(a_start, a_end - 1, a_start, a_end - 1)
    arr = np.asarray(arr, dtype=float)
    # bin starts in local coords
    n = arr.shape[0]
    bin_starts = np.arange(n) * resolution + start
    return arr, bin_starts


def render_pyramid(ax, matrix, bin_starts, resolution, axis_min, axis_max,
                   flip=False, vmax=None, height_frac=0.5):
    """Render an upper-triangle pyramid heatmap via pcolormesh on a rotated
    coordinate grid. The diagonal runs along the x axis. height_frac caps
    how far up the pyramid extends as a fraction of (axis_max-axis_min)."""
    if matrix.size == 0:
        return
    n = matrix.shape[0]
    # Replace zeros with NaN so they render as transparent / blank
    m = matrix.astype(float).copy()
    # symmetrize and keep upper triangle (avoid double-count when mirrored)
    m = np.where(np.isfinite(m), m, 0.0)
    m = (m + m.T) * 0.5
    # mask lower triangle (we only render upper)
    triu = np.triu(np.ones_like(m), k=0).astype(bool)
    plot_m = np.where(triu, m, np.nan)

    # build the rotated mesh of cell corners.
    # cell (i, j) where j >= i has corners at (i, j), (i+1, j), (i+1, j+1), (i, j+1)
    # rotated: corner (a, b) → x = (a + b) * 0.5,  y = (b - a) * 0.5
    a_grid, b_grid = np.meshgrid(np.arange(n + 1), np.arange(n + 1),
                                 indexing="ij")
    xc = (a_grid + b_grid) * 0.5 * resolution + bin_starts[0]
    yc = (b_grid - a_grid) * 0.5 * resolution
    if flip:
        yc = -yc

    finite_vals = plot_m[np.isfinite(plot_m) & (plot_m > 0)]
    if finite_vals.size == 0:
        return
    if vmax is None:
        vmax = np.nanpercentile(finite_vals, 99.0)
    vmin = max(finite_vals.min(), 1.0)
    norm = LogNorm(vmin=vmin, vmax=max(vmax, vmin * 10))
    ax.pcolormesh(xc, yc, plot_m, cmap=WHITE_RED, norm=norm,
                  shading="flat", linewidth=0, antialiased=False,
                  rasterized=True)

    ax.set_xlim(axis_min, axis_max)
    span = axis_max - axis_min
    h = span * height_frac
    if flip:
        ax.set_ylim(-h, 0)
    else:
        ax.set_ylim(0, h)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--ccs", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--hic", required=True)
    ap.add_argument("--hic-chrom", default="assembly")
    ap.add_argument("--hic-offset-a", type=int, required=True)
    ap.add_argument("--hic-offset-b", type=int, required=True)
    ap.add_argument("--hic-resolution", type=int, default=25_000)
    ap.add_argument("--hic-pyramid-frac", type=float, default=None,
                    help="height of each Hi-C pyramid as fraction of "
                         "(window width); smaller = squatter pyramid. "
                         "Default (None) auto-computes the value that makes "
                         "Hi-C cells render as squares for the current panel "
                         "aspect ratio — this is what you want unless you are "
                         "deliberately distorting the heatmap.")
    ap.add_argument("--hapA-scaf", required=True)
    ap.add_argument("--hapA-left", type=int, required=True)
    ap.add_argument("--hapA-right", type=int, required=True)
    ap.add_argument("--hapB-scaf", required=True)
    ap.add_argument("--hapB-left", type=int, required=True)
    ap.add_argument("--hapB-right", type=int, required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--break-window", type=int, default=1000,
                    help="bp; CCS reads whose mapped span overlaps this "
                         "window centered on a breakpoint are highlighted dark")
    ap.add_argument("--break-frac", type=float, default=0.20,
                    help="position of the left breakpoint as a fraction of "
                         "panel width (right breakpoint sits at 1 - this). "
                         "0.20 → breakpoints at 20%% and 80%%.")
    ap.add_argument("--cov-bins", type=int, default=600,
                    help="number of bins for the CCS depth track")
    ap.add_argument("--chain-max-gap", type=int, default=50_000)
    ap.add_argument("--chain-rel-skew", type=float, default=0.3)
    ap.add_argument("--min-input-block", type=int, default=1000)
    ap.add_argument("--min-output-block", type=int, default=5_000)
    ap.add_argument("--title", default=None)
    ap.add_argument("--cartoon-paf", default=None,
                    help="PAF for the top-right chromosome cartoon. If omitted, "
                         "no cartoon is drawn. We use the asm20 PAF here for "
                         "fuller whole-chromosome coverage even though the "
                         "ribbon panel uses asm5.")
    ap.add_argument("--cartoon-min-block", type=int, default=5000)
    ap.add_argument("--zoom-width", type=int, default=40000,
                    help="bp; total width of each breakpoint zoom panel "
                         "(centered on the breakpoint). 40 kb gives ~7.5 kb "
                         "context past a typical 25 kb CCS read.")
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.ccs, "rb")
    fai = {ln.split("\t")[0]: int(ln.split("\t")[1])
           for ln in open(args.fai)}
    hapA_len = fai[args.hapA_scaf]
    hapB_len = fai[args.hapB_scaf]

    # Place the breakpoints at break_frac / (1 - break_frac) of the panel
    # width on both axes. For an inversion span S on either axis, the panel
    # width = S / (1 - 2 * break_frac).
    bf = args.break_frac
    if not 0 < bf < 0.5:
        sys.exit("--break-frac must be in (0, 0.5)")
    ax_lo, ax_hi = sorted([args.hapA_left, args.hapA_right])
    bx_lo, bx_hi = sorted([args.hapB_left, args.hapB_right])
    a_span = ax_hi - ax_lo
    b_span = bx_hi - bx_lo
    a_pad = int(a_span * bf / (1 - 2 * bf))
    b_pad = int(b_span * bf / (1 - 2 * bf))
    xmin = max(0, ax_lo - a_pad); xmax = min(hapA_len, ax_hi + a_pad)
    ymin = max(0, bx_lo - b_pad); ymax = min(hapB_len, bx_hi + b_pad)

    print("  loading + chaining PAF ...", file=sys.stderr)
    chained = load_chained(args.paf, args.hapA_scaf, args.hapB_scaf,
                           args.chain_max_gap, args.chain_rel_skew,
                           args.min_input_block, args.min_output_block)
    in_win = [b for b in chained
              if not (b.t_end < xmin or b.t_start > xmax)
              and not (b.q_end < ymin or b.q_start > ymax)]
    print(f"    {len(in_win)} synteny blocks in window "
          f"(hapA {xmin:,}-{xmax:,} ; hapB {ymin:,}-{ymax:,})",
          file=sys.stderr)

    print("  fetching CCS reads ...", file=sys.stderr)
    A_reads = fetch_reads(bam, args.hapA_scaf, xmin, xmax,
                          (args.hapA_left, args.hapA_right),
                          break_window=args.break_window)
    B_reads = fetch_reads(bam, args.hapB_scaf, ymin, ymax,
                          (args.hapB_left, args.hapB_right),
                          break_window=args.break_window)
    A_cov_x, A_cov_y, _ = fetch_coverage(bam, args.hapA_scaf, xmin, xmax,
                                         args.cov_bins)
    B_cov_x, B_cov_y, _ = fetch_coverage(bam, args.hapB_scaf, ymin, ymax,
                                         args.cov_bins)
    cov_max = max(A_cov_y.max() if A_cov_y.size else 1,
                  B_cov_y.max() if B_cov_y.size else 1, 1.0)

    print("  fetching Hi-C ...", file=sys.stderr)
    hic = hicstraw.HiCFile(args.hic)
    A_mat, A_bins = fetch_hic_pyramid(hic, args.hic_chrom, args.hic_offset_a,
                                      xmin, xmax, args.hic_resolution)
    B_mat, B_bins = fetch_hic_pyramid(hic, args.hic_chrom, args.hic_offset_b,
                                      ymin, ymax, args.hic_resolution)
    print(f"    hapA Hi-C matrix {A_mat.shape}, hapB {B_mat.shape}",
          file=sys.stderr)

    # ------------------- layout -------------------
    # 180 mm wide x 9 in tall figure (single-column MBE width).
    # Anything ≥6 pt remains legible; nothing exceeds the page bounds.
    MM_PER_INCH = 25.4
    fig = plt.figure(figsize=(180 / MM_PER_INCH, 9.0))
    fig.subplots_adjust(left=0.10, right=0.99, top=0.96, bottom=0.05, hspace=0)

    # row heights — Hi-C panels get a reasonable strip, ribbon panel gets
    # more room (it's the most information-dense). Coverage tracks sit
    # outboard of the read piles and are ~half their height.
    rows = [
        ("hicA",   0.30),
        ("covA",   0.022),  # ~75% of previous 0.03
        ("readA",  0.06),
        ("ribbon", 0.13),   # shorter ribbon so chr tracks sit closer to CCS
        ("readB",  0.06),
        ("covB",   0.022),
        ("hicB",   0.30),
    ]
    total = sum(h for _, h in rows)
    # Reserve the bottom ~22% of the figure for the 4-breakpoint zoom row;
    # main panel stack occupies the top ~65% (top edge at y=0.95, hicB bottom
    # ~y=0.30).
    main_height_frac = 0.65
    y_cursor = 0.95
    panels = {}
    for name, h in rows:
        frac = h / total * main_height_frac
        y_cursor -= frac
        panels[name] = fig.add_axes([0.12, y_cursor, 0.86, frac])

    # x-axis range and helpers
    def setup_axis(ax, xlo, xhi):
        ax.set_xlim(xlo, xhi)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

    # --- Hi-C top (hapA) ---
    # Auto-compute the pyramid height fraction so Hi-C cells render as
    # squares: with the pyramid rotated 45° on the (x,y)-plane, the data
    # ranges are span x = (xmax - xmin) and span y = span x * height_frac.
    # For square pixels we need (panel_h_in / panel_w_in) = height_frac.
    bbox = panels["hicA"].get_position()
    panel_w_in = bbox.width  * fig.get_figwidth()
    panel_h_in = bbox.height * fig.get_figheight()
    auto_pyramid_frac = panel_h_in / panel_w_in
    pyramid_frac = (auto_pyramid_frac
                    if args.hic_pyramid_frac is None
                    else args.hic_pyramid_frac)
    render_pyramid(panels["hicA"], A_mat, A_bins, args.hic_resolution,
                   xmin, xmax, flip=False, height_frac=pyramid_frac)
    panels["hicA"].set_title(
        args.title or f"{args.hapA_scaf} × {args.hapB_scaf}",
        fontsize=11, loc="left")
    panels["hicA"].text(-0.005, 0.5, "Hi-C hapA", transform=panels["hicA"].transAxes,
                         ha="right", va="center", fontsize=8, color="#555555")

    # --- hapA coverage (above reads) ---
    cov_ax = panels["covA"]
    setup_axis(cov_ax, xmin, xmax)
    if A_cov_y.size:
        cov_ax.fill_between(A_cov_x, 0, A_cov_y, step="mid",
                            color="#777777", linewidth=0)
    cov_ax.set_ylim(0, cov_max * 1.05)
    for pos in (args.hapA_left, args.hapA_right):
        cov_ax.axvline(pos, color=BREAK_COLOR, ls="--", lw=0.8, alpha=0.7)
    cov_ax.text(-0.005, 0.5, "CCS depth A",
                transform=cov_ax.transAxes,
                ha="right", va="center", fontsize=7, color="#555555")
    cov_ax.text(1.002, 0.95, f"{int(cov_max)}×",
                transform=cov_ax.transAxes,
                ha="left", va="top", fontsize=6, color="#888888")

    # --- hapA reads ---
    setup_axis(panels["readA"], xmin, xmax)
    placed, n_rowsA = stack_reads(A_reads)
    if n_rowsA > 0:
        row_h = 1.0 / n_rowsA
        bg_rects, hi_rects = [], []
        for r, row in placed:
            # leave 15% of the row height as a vertical gap so adjacent
            # stacked reads are visually separate (matters most at small
            # figure sizes where the read pile is dense)
            gap = row_h * 0.15
            ybot = row * row_h + gap / 2
            ytop = (row + 1) * row_h - gap / 2
            rect = [(r["start"], ybot),
                    (r["end"],   ybot),
                    (r["end"],   ytop),
                    (r["start"], ytop)]
            (hi_rects if r["cross"] else bg_rects).append(rect)
        if bg_rects:
            panels["readA"].add_collection(PolyCollection(
                bg_rects, facecolors=READ_COLOR_BG, edgecolors="none",
                linewidths=0, antialiased=False, rasterized=True))
        if hi_rects:
            panels["readA"].add_collection(PolyCollection(
                hi_rects, facecolors=READ_COLOR_CROSS, edgecolors="none",
                linewidths=0, antialiased=False, rasterized=True))
    panels["readA"].set_ylim(0, 1)
    for pos in (args.hapA_left, args.hapA_right):
        panels["readA"].axvline(pos, color=BREAK_COLOR, ls="--", lw=0.8, alpha=0.7)
    panels["readA"].text(-0.005, 0.5, "CCS hapA",
                          transform=panels["readA"].transAxes,
                          ha="right", va="center", fontsize=8, color="#555555")

    # --- ribbon ---
    rb = panels["ribbon"]
    for s in rb.spines.values():
        s.set_visible(False)
    rb.set_xticks([]); rb.set_yticks([])
    rb.set_xlim(0, 1)
    rb.set_ylim(0, 1)

    # Chr tracks sit close to the top/bottom of the (now-shorter) ribbon panel
    # so they hug the read-pile rows just above/below.
    A_Y = 0.82; B_Y = 0.18
    track_h = 0.06
    rb.add_patch(mpatches.Rectangle((0, A_Y - track_h/2), 1, track_h,
                                    facecolor=TRACK_COLOR, edgecolor="none"))
    rb.add_patch(mpatches.Rectangle((0, B_Y - track_h/2), 1, track_h,
                                    facecolor=TRACK_COLOR, edgecolor="none"))

    def xa(p): return (p - xmin) / (xmax - xmin)
    def xb(p): return (p - ymin) / (ymax - ymin)

    fwd_polys = []; rev_polys = []
    for a in in_win:
        ax1, ax2 = xa(a.t_start), xa(a.t_end)
        if a.strand == "+":
            bx1, bx2 = xb(a.q_start), xb(a.q_end)
            fwd_polys.append([(ax1, A_Y - track_h/2),
                              (ax2, A_Y - track_h/2),
                              (bx2, B_Y + track_h/2),
                              (bx1, B_Y + track_h/2)])
        else:
            bx1, bx2 = xb(a.q_end), xb(a.q_start)
            rev_polys.append([(ax1, A_Y - track_h/2),
                              (ax2, A_Y - track_h/2),
                              (bx2, B_Y + track_h/2),
                              (bx1, B_Y + track_h/2)])
    if fwd_polys:
        rb.add_collection(PolyCollection(fwd_polys, facecolors=FWD_COLOR,
                                         edgecolors="none", alpha=0.55))
    if rev_polys:
        rb.add_collection(PolyCollection(rev_polys, facecolors=REV_COLOR,
                                         edgecolors="none", alpha=0.55))

    # Breakpoint markers on ribbon panel + Mb ticks for each track.
    # Each breakpoint also gets a circled-digit label (①②③④) so the four
    # zoom panels at the bottom of the figure can be keyed back to the
    # exact breakpoint they're zooming on.
    CIRCLED = ["①", "②", "③", "④"]  # ① ② ③ ④
    BP_NUMBERS = {
        ("A", "left"):  CIRCLED[0],
        ("A", "right"): CIRCLED[1],
        ("B", "left"):  CIRCLED[2],
        ("B", "right"): CIRCLED[3],
    }
    for side, pos in (("left", args.hapA_left), ("right", args.hapA_right)):
        rb.plot([xa(pos), xa(pos)], [A_Y - track_h/2 - 0.05, A_Y + track_h/2 + 0.05],
                color=BREAK_COLOR, lw=1, ls="--", alpha=0.8, zorder=5)
        rb.text(xa(pos), A_Y + track_h/2 + 0.11, BP_NUMBERS[("A", side)],
                ha="center", va="center", fontsize=12, color="#222222",
                fontweight="bold", zorder=20,
                bbox=dict(boxstyle="circle,pad=0.05",
                          facecolor="white", edgecolor="none", alpha=0.85))
    for side, pos in (("left", args.hapB_left), ("right", args.hapB_right)):
        rb.plot([xb(pos), xb(pos)], [B_Y - track_h/2 - 0.05, B_Y + track_h/2 + 0.05],
                color=BREAK_COLOR, lw=1, ls="--", alpha=0.8, zorder=5)
        rb.text(xb(pos), B_Y - track_h/2 - 0.11, BP_NUMBERS[("B", side)],
                ha="center", va="center", fontsize=12, color="#222222",
                fontweight="bold", zorder=20,
                bbox=dict(boxstyle="circle,pad=0.05",
                          facecolor="white", edgecolor="none", alpha=0.85))

    def add_ticks(y, lo, hi, fn, above):
        vals = list(np.linspace(lo, hi, 7))
        for i, v in enumerate(vals):
            rb.plot([fn(v), fn(v)],
                    [y + 0.015, y + 0.04] if above else [y - 0.015, y - 0.04],
                    color=TRACK_COLOR, lw=0.6)
            # Mark the rightmost tick with the "Mbp" unit so the units of
            # the chr-track axis are unambiguous.
            txt = f"{v/1e6:.3f}"
            if i == len(vals) - 1:
                txt += " Mbp"
            rb.text(fn(v), y + 0.055 if above else y - 0.055,
                    txt, ha="center",
                    va="bottom" if above else "top",
                    fontsize=6, color=TRACK_COLOR)
    add_ticks(A_Y, xmin, xmax, xa, above=True)
    add_ticks(B_Y, ymin, ymax, xb, above=False)
    def short_scaf(name):
        if name.startswith("eupHapAv0.3_"):
            return "Hap A " + name.replace("eupHapAv0.3_", "")
        if name.startswith("eupHapBv0.3_"):
            return "Hap B " + name.replace("eupHapBv0.3_", "")
        return name
    rb.text(-0.005, A_Y, short_scaf(args.hapA_scaf), transform=rb.transAxes,
            ha="right", va="center", fontsize=8)
    rb.text(-0.005, B_Y, short_scaf(args.hapB_scaf), transform=rb.transAxes,
            ha="right", va="center", fontsize=8)
    # Legend goes in the strip *just above the 4 zoom panels*, right-justified.
    # zoom_y_top is ~0.22, so anchoring the legend at y=0.295 (above the zoom
    # panel titles which extend slightly above the panel boxes) keeps it from
    # overlapping the per-panel titles.
    fig.legend(handles=[
        mpatches.Patch(color=FWD_COLOR, label="alignment (+) colinear"),
        mpatches.Patch(color=REV_COLOR, label="alignment (−) inverted"),
        mpatches.Patch(color=READ_COLOR_CROSS, label="CCS read spans breakpoint"),
        mpatches.Patch(color=BREAK_COLOR, label="breakpoint"),
    ], fontsize=6.5, loc="lower right", bbox_to_anchor=(0.99, 0.275),
       ncol=2, frameon=True, framealpha=0.92, edgecolor="#bbbbbb",
       handlelength=1.5, columnspacing=0.8, handletextpad=0.4)

    # --- inversion-locator cartoon (top-right, above the legend) ---
    if args.cartoon_paf:
        cartoon_ax = fig.add_axes([0.62, 0.86, 0.36, 0.085])
        draw_chr_cartoon(cartoon_ax, args.cartoon_paf,
                         args.hapA_scaf, args.hapB_scaf, hapA_len, hapB_len,
                         args.hapA_left, args.hapA_right,
                         args.hapB_left, args.hapB_right,
                         min_block=args.cartoon_min_block, chain=True)

    # --- hapB reads ---
    setup_axis(panels["readB"], ymin, ymax)
    placed, n_rowsB = stack_reads(B_reads)
    if n_rowsB > 0:
        row_h = 1.0 / n_rowsB
        bg_rects, hi_rects = [], []
        for r, row in placed:
            gap = row_h * 0.15
            ytop = 1 - row * row_h - gap / 2
            ybot = 1 - (row + 1) * row_h + gap / 2
            rect = [(r["start"], ybot),
                    (r["end"],   ybot),
                    (r["end"],   ytop),
                    (r["start"], ytop)]
            (hi_rects if r["cross"] else bg_rects).append(rect)
        if bg_rects:
            panels["readB"].add_collection(PolyCollection(
                bg_rects, facecolors=READ_COLOR_BG, edgecolors="none",
                linewidths=0, antialiased=False, rasterized=True))
        if hi_rects:
            panels["readB"].add_collection(PolyCollection(
                hi_rects, facecolors=READ_COLOR_CROSS, edgecolors="none",
                linewidths=0, antialiased=False, rasterized=True))
    panels["readB"].set_ylim(0, 1)
    for pos in (args.hapB_left, args.hapB_right):
        panels["readB"].axvline(pos, color=BREAK_COLOR, ls="--", lw=0.8, alpha=0.7)
    panels["readB"].text(-0.005, 0.5, "CCS hapB",
                          transform=panels["readB"].transAxes,
                          ha="right", va="center", fontsize=8, color="#555555")

    # --- hapB coverage (below reads, mirrored) ---
    cov_ax = panels["covB"]
    setup_axis(cov_ax, ymin, ymax)
    if B_cov_y.size:
        cov_ax.fill_between(B_cov_x, 0, B_cov_y, step="mid",
                            color="#777777", linewidth=0)
    cov_ax.set_ylim(cov_max * 1.05, 0)  # inverted so depth points down
    for pos in (args.hapB_left, args.hapB_right):
        cov_ax.axvline(pos, color=BREAK_COLOR, ls="--", lw=0.8, alpha=0.7)
    cov_ax.text(-0.005, 0.5, "CCS depth B",
                transform=cov_ax.transAxes,
                ha="right", va="center", fontsize=7, color="#555555")
    cov_ax.text(1.002, 0.05, f"{int(cov_max)}×",
                transform=cov_ax.transAxes,
                ha="left", va="bottom", fontsize=6, color="#888888")

    # --- Hi-C bottom (hapB), upside-down ---
    render_pyramid(panels["hicB"], B_mat, B_bins, args.hic_resolution,
                   ymin, ymax, flip=True, height_frac=pyramid_frac)
    panels["hicB"].text(-0.005, 0.5, "Hi-C hapB",
                         transform=panels["hicB"].transAxes,
                         ha="right", va="center", fontsize=8, color="#555555")

    # --- breakpoint zoom row (bottom of figure) ---
    zw = args.zoom_width
    zw_half = zw // 2
    # Each breakpoint definition: circled label, panel suffix, scaffold, bp,
    # the haplotype's max length (for clipping), and the title colour.
    bp_specs = [
        (CIRCLED[0], "hapA left",  args.hapA_scaf, args.hapA_left,  hapA_len, "#1f4f8a"),
        (CIRCLED[1], "hapA right", args.hapA_scaf, args.hapA_right, hapA_len, "#1f4f8a"),
        (CIRCLED[2], "hapB left",  args.hapB_scaf, args.hapB_left,  hapB_len, "#226a1f"),
        (CIRCLED[3], "hapB right", args.hapB_scaf, args.hapB_right, hapB_len, "#226a1f"),
    ]
    # 4 panels horizontally, moved up just under the main panel stack
    # (main panels currently bottom-out at y ~ 0.25 with main_height_frac=0.70)
    zoom_y_top = 0.22
    zoom_y_bot = 0.05
    zoom_h    = zoom_y_top - zoom_y_bot
    zoom_left = 0.05
    zoom_right = 0.99
    zoom_gap_frac = 0.030      # bigger horizontal gap between zoom panels
    panel_w = ((zoom_right - zoom_left) - 3 * zoom_gap_frac) / 4

    # Section header on the LEFT (legend sits on the right at the same height).
    fig.text(zoom_left, zoom_y_top + 0.060,
             "Breakpoint zooms (CCS reads, ±{:.0f} kb)".format(zw / 2000),
             ha="left", va="bottom", fontsize=9, fontweight="bold",
             color="#333333")

    for i, (circle, label, scaf, bp, scaf_len, title_color) in enumerate(bp_specs):
        x0 = zoom_left + i * (panel_w + zoom_gap_frac)
        zoom_ax = fig.add_axes([x0, zoom_y_bot, panel_w, zoom_h])
        z_lo = max(0, bp - zw_half)
        z_hi = min(scaf_len, bp + zw_half)
        z_reads = fetch_reads(bam, scaf, z_lo, z_hi, (bp,),
                              break_window=args.break_window)
        placed, n_rows_z = stack_reads(z_reads)
        if n_rows_z > 0:
            row_h = 1.0 / n_rows_z
            bg_rects, hi_rects = [], []
            for r, row in placed:
                gap = row_h * 0.15
                ybot = row * row_h + gap / 2
                ytop = (row + 1) * row_h - gap / 2
                rect = [(r["start"], ybot),
                        (r["end"],   ybot),
                        (r["end"],   ytop),
                        (r["start"], ytop)]
                (hi_rects if r["cross"] else bg_rects).append(rect)
            if bg_rects:
                zoom_ax.add_collection(PolyCollection(
                    bg_rects, facecolors=READ_COLOR_BG, edgecolors="none",
                    linewidths=0, antialiased=False, rasterized=True))
            if hi_rects:
                zoom_ax.add_collection(PolyCollection(
                    hi_rects, facecolors=READ_COLOR_CROSS, edgecolors="none",
                    linewidths=0, antialiased=False, rasterized=True))
        zoom_ax.axvline(bp, color=BREAK_COLOR, ls="--", lw=1.0, alpha=0.85)
        zoom_ax.set_xlim(z_lo, z_hi)
        zoom_ax.set_ylim(0, 1)
        zoom_ax.set_yticks([])
        zoom_ax.tick_params(axis="x", labelsize=6, length=2, pad=1)
        # Two ticks (panel edges) — the breakpoint position itself is
        # already given in the panel title and as the dashed red line, so
        # the middle tick would be redundant. Tick labels are rotated 30°
        # so the rightmost tick of one panel can't horizontally overlap
        # the leftmost tick of the next panel at this 180 mm width.
        ticks = np.linspace(z_lo, z_hi, 2)
        zoom_ax.set_xticks(ticks)
        zoom_ax.set_xticklabels([f"{t/1e6:.3f}" for t in ticks],
                                rotation=30, ha="right",
                                rotation_mode="anchor")
        zoom_ax.tick_params(axis="x", which="major", pad=1)
        zoom_ax.set_xlabel("Mbp", fontsize=6.5, labelpad=1)
        for s in ("top", "right", "left"):
            zoom_ax.spines[s].set_visible(False)
        zoom_ax.spines["bottom"].set_color("#888888")
        # 1 kb scale bar in the top-right of the zoom panel
        scale_x1 = z_hi - 0.1 * zw
        scale_x0 = scale_x1 - 1000
        zoom_ax.plot([scale_x0, scale_x1], [0.96, 0.96],
                     color="#000000", lw=1.2, solid_capstyle="butt")
        zoom_ax.text((scale_x0 + scale_x1) / 2, 0.99, "1 kb",
                     ha="center", va="bottom", fontsize=6, color="#333333")
        # Two-line title: big circled-number key on the left, panel
        # descriptor on the right (broken across two lines so it fits in
        # the narrow 180/4 mm panel width without truncation).
        scaf_short = scaf.replace('eupHapAv0.3_', '').replace('eupHapBv0.3_', '')
        zoom_ax.text(-0.005, 1.05, circle,
                     transform=zoom_ax.transAxes,
                     ha="left", va="bottom",
                     fontsize=13, fontweight="bold", color=title_color)
        zoom_ax.text(0.11, 1.27,
                     f"{label} ({scaf_short})",
                     transform=zoom_ax.transAxes,
                     ha="left", va="top", fontsize=7,
                     color=title_color, fontweight="bold")
        zoom_ax.text(0.11, 1.07,
                     f"{bp/1e6:.3f} Mbp · n={len(z_reads)} reads",
                     transform=zoom_ax.transAxes,
                     ha="left", va="top", fontsize=6.5, color="#444444")

    out_pdf = args.out_prefix + ".style4.pdf"
    out_png = args.out_prefix + ".style4.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=180)
    print(f"wrote {out_pdf}", file=sys.stderr)
    print(f"wrote {out_png}", file=sys.stderr)


if __name__ == "__main__":
    main()

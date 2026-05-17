#!/usr/bin/env python3
"""Driver: render one ``hapduo-style4`` panel per row of a breakpoints TSV.

The TSV must have the same schema as the output of ``hapduo-detect``:

    rank  hapA_chr  hapA_pos  hapB_chr  hapB_pos  end  inv_len  n_anchors

with two rows per inversion (``end=left`` and ``end=right``). Use the same
script repeatedly: pair the two rows per rank and call
``hapduo-style4`` once per inversion with the right command-line flags.

Hi-C contact-map offsets within a competitively-mapped ``.hic`` file are
computed automatically from the FASTA index (``.fai``).
"""
from __future__ import annotations
import argparse
import csv
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def ordinal(n):
    if 10 <= n % 100 <= 20:
        return f"{n}th"
    return f"{n}{ {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th') }"


def load_fai_offsets(fai_path):
    offsets = {}
    running = 0
    for ln in open(fai_path):
        name, length = ln.split()[0], int(ln.split()[1])
        offsets[name] = running
        running += length
    return offsets


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True,
                    help="breakpoints TSV (hapduo-detect output)")
    ap.add_argument("--paf", required=True,
                    help="asm5 hapB-to-hapA PAF (used for ribbon panel + read panels)")
    ap.add_argument("--cartoon-paf", default=None,
                    help="asm20 PAF for the chromosome-scale cartoon (defaults to --paf)")
    ap.add_argument("--ccs", default=None,
                    help="coordinate-sorted CCS BAM (optional; required for read pile + depth)")
    ap.add_argument("--fai", required=True,
                    help="FASTA index of the concatenated assembly (used for chr offsets)")
    ap.add_argument("--hic", default=None,
                    help="Juicer .hic file (optional; required for Hi-C panels)")
    ap.add_argument("--hic-chrom", default="assembly",
                    help='chromosome name inside the .hic (default "assembly" for '
                         'competitively-mapped single-pseudo-contig setups)')
    ap.add_argument("--hic-resolution", type=int, default=25000)
    ap.add_argument("--outdir", required=True,
                    help="output directory for rendered figures")
    ap.add_argument("--top", type=int, default=None,
                    help="render only the top N inversions by inv_len (default: all)")
    ap.add_argument("--style4-script", default=None,
                    help="path to plot_style4_hic.py; default uses hapduo-style4 from PATH")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fai_offsets = load_fai_offsets(args.fai)
    cartoon_paf = args.cartoon_paf or args.paf

    # Decide how to invoke the per-inversion renderer
    if args.style4_script:
        invoker = ["python3", args.style4_script]
    elif shutil.which("hapduo-style4"):
        invoker = ["hapduo-style4"]
    else:
        invoker = ["python3", "-m", "hapduo.plot_style4_hic"]

    by_rank = defaultdict(dict)
    with open(args.tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_rank[int(row["rank"])][row["end"]] = row

    ranks = sorted(by_rank.keys())
    if args.top:
        ranks = ranks[:args.top]
    print(f"# rendering {len(ranks)} inversions: ranks {ranks}", file=sys.stderr)

    for rank in ranks:
        entry = by_rank[rank]
        if "left" not in entry or "right" not in entry:
            print(f"!! rank {rank} is missing left/right row; skipping", file=sys.stderr)
            continue
        L, R = entry["left"], entry["right"]
        hapA_scaf = L["hapA_chr"]; hapB_scaf = L["hapB_chr"]
        hapA_left  = int(L["hapA_pos"]); hapA_right = int(R["hapA_pos"])
        hapB_left  = int(L["hapB_pos"]); hapB_right = int(R["hapB_pos"])
        inv_len = int(L["inv_len"]); n_anch = int(L["n_anchors"])

        chr_token = hapA_scaf.split("_")[-1]
        prefix = outdir / f"rank{rank:02d}_{chr_token}"
        chr_num = chr_token.replace("chr", "")
        title = (f"{ordinal(rank)} largest in genome — chr {chr_num}  "
                 f"({inv_len/1e6:.2f} Mb, {n_anch} anchors)")

        cmd = invoker + [
            "--paf", args.paf,
            "--fai", args.fai,
            "--hapA-scaf", hapA_scaf,
            "--hapA-left", str(hapA_left),
            "--hapA-right", str(hapA_right),
            "--hapB-scaf", hapB_scaf,
            "--hapB-left", str(hapB_left),
            "--hapB-right", str(hapB_right),
            "--cartoon-paf", cartoon_paf,
            "--out-prefix", str(prefix),
            "--title", title,
        ]
        if args.ccs:
            cmd += ["--ccs", args.ccs]
        if args.hic and hapA_scaf in fai_offsets and hapB_scaf in fai_offsets:
            cmd += ["--hic", args.hic,
                    "--hic-chrom", args.hic_chrom,
                    "--hic-offset-a", str(fai_offsets[hapA_scaf]),
                    "--hic-offset-b", str(fai_offsets[hapB_scaf]),
                    "--hic-resolution", str(args.hic_resolution)]

        print(f"\n# rank {rank} -> {prefix}.style4.{{pdf,png}}", file=sys.stderr)
        r = subprocess.run(cmd)
        if r.returncode != 0:
            print(f"!! rank {rank} failed (exit {r.returncode})", file=sys.stderr)


if __name__ == "__main__":
    main()

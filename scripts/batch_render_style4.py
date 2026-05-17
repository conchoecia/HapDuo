#!/usr/bin/env python3
"""Driver: render Style-4 plots for all rows of breakpoints_v4.tsv.

For each rank, pair the 'left' and 'right' rows to get
(hapA_scaf, hapA_left, hapA_right, hapB_scaf, hapB_left, hapB_right),
then invoke plot_style4_hic.py.

Run from /media/darrin/euplokamis_inv_work.
"""
from __future__ import annotations
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

H03  = Path(__file__).resolve().parent.parent          # hypotheses/H03_.../
TSV  = H03 / "breakpoints_v5.tsv"
SCRIPT = H03 / "scripts" / "plot_style4_hic.py"
OUTDIR = H03 / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

PAF        = str(H03 / "work" / "hapB_to_hapA.paf")
CARTOON_PAF= str(H03 / "work" / "hapB_to_hapA.asm20.paf")
CCS        = str(H03 / "work" / "ccs_to_hapAhapB.sorted.bam")
FAI        = "/media/darrin/euplokamis_reads_data/new_folder/v0.3/eup_hapAhapBMetaExtraMito_v0.3.fasta.fai"
HIC        = "/media/darrin/euplokamis_reads_data/new_folder/v0.3/hic/GAP_hic_map5/output/euplov03/euplov03_allLibs_q_5.hic"
HIC_CHROM  = "assembly"

# Pre-computed bp offsets within the concatenated 'assembly' pseudo-contig.
fai_offsets = {}
running = 0
for ln in open(FAI):
    name, length = ln.split()[0], int(ln.split()[1])
    fai_offsets[name] = running
    running += length

# Parse breakpoints_v4.tsv -> { rank: {'left': row, 'right': row} }
by_rank = defaultdict(dict)
with open(TSV) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        by_rank[int(row["rank"])][row["end"]] = row

ranks = sorted(by_rank.keys())[:10]
print(f"# rendering ranks {ranks}", file=sys.stderr)

for rank in ranks:
    entry = by_rank[rank]
    left  = entry["left"]
    right = entry["right"]
    hapA_scaf = left["hapA_chr"]
    hapB_scaf = left["hapB_chr"]
    hapA_left  = int(left["hapA_pos"])
    hapA_right = int(right["hapA_pos"])
    hapB_left  = int(left["hapB_pos"])
    hapB_right = int(right["hapB_pos"])
    inv_len    = int(left["inv_len"])
    n_anch     = int(left["n_anchors"])

    chr_token = hapA_scaf.split("_")[-1]
    prefix = OUTDIR / f"rank{rank:02d}_{chr_token}"

    def ordinal(n):
        if 10 <= n % 100 <= 20:
            return f"{n}th"
        return f"{n}{ {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th') }"

    chr_num = chr_token.replace("chr", "")
    title = (f"{ordinal(rank)} largest in genome — chr {chr_num}  "
             f"({inv_len/1e6:.2f} Mb, {n_anch} anchors)")

    cmd = [
        "python3", str(SCRIPT),
        "--paf", PAF,
        "--ccs", CCS,
        "--fai", FAI,
        "--hic", HIC,
        "--hic-chrom", HIC_CHROM,
        "--hic-offset-a", str(fai_offsets[hapA_scaf]),
        "--hic-offset-b", str(fai_offsets[hapB_scaf]),
        "--hic-resolution", "25000",
        "--hapA-scaf", hapA_scaf,
        "--hapA-left", str(hapA_left),
        "--hapA-right", str(hapA_right),
        "--hapB-scaf", hapB_scaf,
        "--hapB-left", str(hapB_left),
        "--hapB-right", str(hapB_right),
        "--cartoon-paf", CARTOON_PAF,
        "--out-prefix", str(prefix),
        "--title", title,
    ]
    print(f"\n# rank {rank} -> {prefix}.style4.{{pdf,png}}", file=sys.stderr)
    print("# " + " ".join(cmd), file=sys.stderr)
    r = subprocess.run(cmd)
    if r.returncode != 0:
        print(f"!! rank {rank} failed (exit {r.returncode})", file=sys.stderr)

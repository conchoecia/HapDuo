# DPGB inter-haplotype inversion plotting scripts

Six Python scripts contributed from the *Euplokamis dunlapae* genome project
(Schultz et al., in prep.) for plotting inter-haplotype inversion calls
from a base-level PAF.

## Pipeline overview

```
hapB.fasta + hapA.fasta
        │
        ▼  minimap2 -cx asm5  (or asm20)
   hapB_to_hapA.paf
        │
        ▼  _paf_chain.py       (greedy chainer)
   chained syntenic blocks
        │
        ├──▶ find_inversions_strandswitch.py
        │      detects inversions by the canonical
        │      +++ −−− +++ strand-switch pattern with
        │      antiparallel-chain trimming, intrusion
        │      tolerance, and Jaccard dedup.
        │      OUTPUT: breakpoints.tsv
        │
        ├──▶ plot_synteny_ideogram.py
        │      whole-chromosome two-haplotype synteny ribbon
        │      ideogram, one row per chromosome pair
        │
        ├──▶ plot_inversion_ideogram.py
        │      ideogram of inversion calls along the haplotype-A
        │      chromosomes, color-coded by size
        │
        └──▶ plot_style4_hic.py
               full multi-evidence inversion-validation panel
               (Hi-C + CCS depth + read pile + ribbon + chr cartoon
               + 4-breakpoint zoom row)

batch_render_style4.py is a driver that reads breakpoints.tsv and
invokes plot_style4_hic.py per inversion.
```

## Key script: `plot_style4_hic.py`

Produces a single 180 × 228 mm figure with:

1. **Top:** Hi-C contact pyramid heatmap for the haplotype-A window
   (white-to-red log scale, auto-computed cell aspect so cells render
   square at any figure size).
2. CCS read-depth track + read pile on haplotype-A, with reads spanning
   a 1 kb window around any breakpoint coloured dark grey.
3. Haplotype-A chromosome track with dashed red breakpoint markers
   labelled ① (left) / ② (right).
4. **Centre:** ribbon panel rendering each chained PAF block as a
   coloured polygon (blue `+` colinear, red `−` inverted), with strand
   transitions placed at fixed 20% and 80% of panel width.
5. Haplotype-B chromosome track with breakpoint markers ③ / ④.
6. CCS read pile + depth on haplotype-B.
7. **Bottom:** mirrored Hi-C pyramid for haplotype-B (y-axis inverted).
8. **Top-right inset:** chromosome-scale two-haplotype cartoon
   (asm20 PAF) showing both haplotypes at true Mb length with this
   inversion's location outlined as a black rectangle on both bars.
9. **Bottom row:** four ±20 kb (40 kb total) CCS read-pile zooms,
   one per breakpoint, keyed by the same circled digits ①②③④ used on
   the chromosome tracks above.

All text is ≥ 6 pt and the PDF is saved with `pdf.fonttype = 42` so
text remains TrueType-vector (editable in Illustrator).

## Other scripts

- `plot_synteny_ideogram.py` — whole-genome two-haplotype synteny
  ideogram (one row per chromosome pair; blue + colinear / red −
  inverted ribbons).
- `plot_inversion_ideogram.py` — single-haplotype ideogram of called
  inversions only, colour-coded by size.
- `_paf_chain.py` — greedy same-strand PAF chainer (target/query gap
  ≤ max_gap, relative skew ≤ rel_skew).
- `find_inversions_strandswitch.py` — antiparallel-chain inversion
  detector with `+` intrusion tolerance and Jaccard dedup.
- `batch_render_style4.py` — driver: reads a breakpoints TSV with
  `rank/end/hapA_chr/hapA_pos/hapB_chr/hapB_pos/inv_len/n_anchors`
  columns and renders one Style 4 figure per rank.

## Dependencies

```
python >= 3.9
numpy
matplotlib >= 3.9
pysam               # for CCS BAM
hic-straw >= 1.3    # for .hic Hi-C contact extraction
```

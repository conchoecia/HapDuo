# HapDuo

**Multi-evidence visualisation of structural variants between two haplotype assemblies of the same individual.**

HapDuo aligns two phased haplotype assemblies, detects inversions from the
canonical `+++ −−− +++` strand-switch pattern in a chained PAF, and renders
a stack of supporting figures — including a per-inversion multi-evidence
panel that combines Hi-C contact pyramids, CCS read piles, PAF synteny
ribbons, a chromosome-scale cartoon, and close-up zoom panels around all
four breakpoints, all on a single page.

It was developed alongside the *Euplokamis dunlapae* genome paper (Schultz
et al., in prep.) where it produced Supplementary Figures S10–S19.

## Install

```bash
pip install hapduo
```

Requires Python ≥ 3.9, `numpy`, `matplotlib`, `pysam`, and `hic-straw`. The
Snakemake pipeline additionally needs `minimap2` and `snakemake` on `PATH`.

## What you get

After install you have five command-line tools:

| Command           | Output |
|-------------------|--------|
| `hapduo-detect`   | Inversion breakpoint TSV from a chained PAF |
| `hapduo-style4`   | One multi-evidence per-inversion panel (180 × 228 mm, ≥ 6 pt text, TrueType — Illustrator-editable) |
| `hapduo-batch`    | All `hapduo-style4` panels for one breakpoints TSV |
| `hapduo-ideogram` | Single-haplotype ideogram of inversion calls |
| `hapduo-synteny`  | Whole-genome two-haplotype synteny ribbon ideogram |

…plus a top-level `Snakefile` that runs the whole pipeline from a single
config file.

## Quickstart (manual)

```bash
# 1. Align hapB → hapA at two sensitivity levels
minimap2 -cx asm5  --eqx -t 16 hapA.fasta hapB.fasta > hapB_to_hapA.asm5.paf
minimap2 -cx asm20 --eqx -t 16 hapA.fasta hapB.fasta > hapB_to_hapA.asm20.paf

# 2. Build a FASTA index of a concatenated hapA + hapB reference (for chr offsets)
cat hapA.fasta hapB.fasta > both.fasta && samtools faidx both.fasta

# 3. Call inversions
hapduo-detect hapB_to_hapA.asm5.paf --min-anchor 5000 --min-inv 50000 > breakpoints.tsv

# 4. Render the whole-genome synteny ideogram (asm20 gives fuller coverage)
hapduo-synteny --paf hapB_to_hapA.asm20.paf --fai both.fasta.fai \
               --chain --min-block 5000 --out-prefix figures/synteny

# 5. Render one multi-evidence panel per inversion
hapduo-batch --tsv breakpoints.tsv \
             --paf hapB_to_hapA.asm5.paf --cartoon-paf hapB_to_hapA.asm20.paf \
             --ccs ccs.sorted.bam --hic contacts.hic --hic-chrom assembly \
             --fai both.fasta.fai --outdir figures/
```

## Quickstart (Snakemake)

```bash
git clone https://github.com/conchoecia/HapDuo.git
cd HapDuo
cp config.example.yaml config.yaml          # edit paths
snakemake --cores 16
```

The pipeline runs steps 1–5 above and writes everything under `outdir/`.

## What does each figure show?

* **Synteny ideogram** (`hapduo-synteny`): one row per chromosome pair,
  haplotype-A bar on top, haplotype-B bar on bottom (both at true Mb
  length), every chained PAF block overlaid as a polygon ribbon
  (blue `+` colinear, red `−` inverted). Quickest way to see where the
  inversions are along the genome.

* **Inversion ideogram** (`hapduo-ideogram`): single-haplotype chromosome
  bars with each inversion call drawn as a coloured rectangle, sized by
  inversion length. Counts per chromosome are summarised in a side panel.

* **Per-inversion multi-evidence panel** (`hapduo-style4`): for a single
  inversion of interest, a single 180 × 228 mm panel showing
    1. Hi-C contact pyramid for the haplotype-A window
       (cells auto-aspected so they render as squares at any figure size).
    2. CCS read-depth track and read pile on haplotype-A, with reads
       spanning a 1 kb window around any breakpoint highlighted dark.
    3. Haplotype-A chromosome track with the two breakpoints marked
       and keyed by ① and ②.
    4. The ribbon panel itself with chained PAF polygons.
    5. Haplotype-B chromosome track with markers ③ and ④.
    6. CCS read pile + depth on haplotype-B.
    7. Mirrored Hi-C pyramid for haplotype-B.
    8. *Top-right inset:* chromosome-scale two-haplotype synteny cartoon
       with the inversion outlined as a black rectangle on both bars.
    9. *Bottom row:* four ±20 kb CCS read-pile zooms, one per breakpoint,
       keyed by the same ①②③④ digits used above.

## Citing

If you use HapDuo, please cite the *Euplokamis dunlapae* genome paper
(Schultz et al., *in prep.*) — the citation will be updated here when the
paper is out.

## Releasing

HapDuo is published to PyPI via OIDC trusted publishing — no API tokens
are stored in the repository. To cut a release:

```bash
git tag v0.1.0
git push --tags
```

The `Release to PyPI` workflow at `.github/workflows/release.yml` builds
an sdist + wheel from the tagged commit, runs a smoke test against the
five console scripts, publishes to PyPI, and attaches the dist files to
a matching GitHub release.

The one-time setup on the PyPI side is documented at the top of that
workflow file.

## Legacy: DPGB dotplot pipeline

This repository was previously called **DPGB** (Dot Plot Genome Browser)
and shipped a Snakemake pipeline for chained-PAF dot plots of CLR + CCS
read mappings. That pipeline is preserved under
[`archive/dpgb-dotplot/`](archive/dpgb-dotplot/) for anyone still running
it; it is not maintained.

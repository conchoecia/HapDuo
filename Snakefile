# HapDuo pipeline: align two haplotype assemblies, detect inter-haplotype
# inversions, and render the supporting multi-evidence figures.
#
# Inputs (declared in config.yaml):
#   hapA_fasta, hapB_fasta              # required: the two assemblies
#   fai                                  # required: FASTA index of a
#                                        #   concatenated reference (hapA +
#                                        #   hapB; built once by the user)
#   ccs_bam                              # optional: HiFi reads aligned to
#                                        #   the concatenated reference
#   hic                                  # optional: .hic contact map of the
#                                        #   concatenated reference
#   hic_chrom (default "assembly")       # optional: chromosome name inside
#                                        #   the .hic for competitive mapping
#   outdir                               # required: output directory
#   threads (default 16)                 # optional: minimap2 thread count
#
# Outputs:
#   {outdir}/paf/hapB_to_hapA.asm5.paf
#   {outdir}/paf/hapB_to_hapA.asm20.paf
#   {outdir}/breakpoints.tsv
#   {outdir}/figures/synteny_ideogram.{pdf,png}
#   {outdir}/figures/inversion_ideogram.{pdf,png}
#   {outdir}/figures/rank{NN}_{chr}.style4.{pdf,png}   # one per inversion
#
# Usage:
#   snakemake --configfile config.yaml --cores 16
#
configfile: "config.yaml"

OUTDIR     = config["outdir"]
HAPA       = config["hapA_fasta"]
HAPB       = config["hapB_fasta"]
FAI        = config["fai"]
THREADS    = config.get("threads", 16)
CCS        = config.get("ccs_bam")
HIC        = config.get("hic")
HIC_CHROM  = config.get("hic_chrom", "assembly")
TOP        = config.get("top", None)

PAF5  = f"{OUTDIR}/paf/hapB_to_hapA.asm5.paf"
PAF20 = f"{OUTDIR}/paf/hapB_to_hapA.asm20.paf"
TSV   = f"{OUTDIR}/breakpoints.tsv"

rule all:
    input:
        f"{OUTDIR}/figures/synteny_ideogram.pdf",
        f"{OUTDIR}/figures/inversion_ideogram.pdf",
        f"{OUTDIR}/figures/.batch_done"

rule minimap2_asm5:
    input:  hapA=HAPA, hapB=HAPB
    output: PAF5
    threads: THREADS
    shell:
        "minimap2 -cx asm5 --eqx -t {threads} {input.hapA} {input.hapB} > {output}"

rule minimap2_asm20:
    input:  hapA=HAPA, hapB=HAPB
    output: PAF20
    threads: THREADS
    shell:
        "minimap2 -cx asm20 --eqx -t {threads} {input.hapA} {input.hapB} > {output}"

rule detect_inversions:
    input:  PAF5
    output: TSV
    shell:
        "hapduo-detect {input} --min-anchor 5000 --min-inv 50000 --top 100 > {output}"

rule synteny_ideogram:
    input:  paf=PAF20, fai=FAI
    output: pdf = f"{OUTDIR}/figures/synteny_ideogram.pdf",
            png = f"{OUTDIR}/figures/synteny_ideogram.png"
    params: prefix = f"{OUTDIR}/figures/synteny_ideogram"
    shell:
        "hapduo-synteny --paf {input.paf} --fai {input.fai} --chain "
        "--min-block 5000 --out-prefix {params.prefix}"

rule inversion_ideogram:
    input:  tsv=TSV, fai=FAI
    output: pdf = f"{OUTDIR}/figures/inversion_ideogram.pdf",
            png = f"{OUTDIR}/figures/inversion_ideogram.png"
    params: prefix = f"{OUTDIR}/figures/inversion_ideogram"
    shell:
        "hapduo-ideogram --breakpoints {input.tsv} --fai {input.fai} "
        "--min-anchors 3 --out-prefix {params.prefix}"

rule batch_style4:
    """Render one multi-evidence per-inversion panel per row of breakpoints.tsv.

    Hi-C and CCS are passed through if available in the config, otherwise the
    per-inversion script renders the panels it can.
    """
    input:  tsv=TSV, paf5=PAF5, paf20=PAF20, fai=FAI
    output: touch(f"{OUTDIR}/figures/.batch_done")
    params:
        outdir = f"{OUTDIR}/figures",
        ccs    = f"--ccs {CCS}" if CCS else "",
        hic    = (f"--hic {HIC} --hic-chrom {HIC_CHROM}" if HIC else ""),
        top    = (f"--top {TOP}" if TOP else "")
    shell:
        "hapduo-batch --tsv {input.tsv} --paf {input.paf5} "
        "--cartoon-paf {input.paf20} --fai {input.fai} "
        "{params.ccs} {params.hic} --outdir {params.outdir} {params.top}"

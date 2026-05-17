"""HapDuo — multi-evidence visualisation of structural variants between
two haplotype assemblies of the same individual.

The package exposes five scripts as console entry points:

    hapduo-detect       chained-PAF strand-switch inversion detector
    hapduo-style4       per-inversion multi-evidence panel
                        (Hi-C + CCS + ribbon + chromosome cartoon + 4 zooms)
    hapduo-ideogram     single-haplotype inversion-call ideogram
    hapduo-synteny      whole-genome two-haplotype synteny ideogram
    hapduo-batch        driver: render one hapduo-style4 figure per row of a
                        breakpoints TSV

…plus a top-level Snakefile that wires the whole pipeline together from a
single config.yaml.
"""

__version__ = "0.1.2"

"""
This gets a list of breakpoints
"""

from Bio import SeqIO
from itertools import islice
from itertools import product
import operator
import pandas as pd
import pickle
import time
import yaml

configfile: "config.yaml"

DNA_comp = {"G": "C",
            "A": "T",
            "C": "G",
            "T": "A"}


"""
https://stackoverflow.com/questions/9475241/split-string-every-nth-character
"""
def split_every(n, iterable):
    i = iter(iterable)
    piece = list(islice(i, n))
    while piece:
        yield piece
        piece = list(islice(i, n))

from operator import itemgetter
from itertools import groupby
def find_ranges(iterable):
    """Yield range of consecutive numbers."""
    ranges = []
    for k,g in groupby(enumerate(iterable),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        ranges.append((group[0],group[-1]))
    return ranges

# how to get the canonical sequence
tab = str.maketrans("ACTG", "TGAC")
def revcomp(seq):
    return seq.translate(tab)[::-1]
def canonical(seq):
    return sorted([seq, revcomp(seq)])[0]

# make all the points
config["points"] = set()
for thisbreak in config["breaks"]:
    thispoint = "-".join([str(x) for x in thisbreak])
    config["points"].add(thispoint)

rule all:
    input:
        expand("synbreaksupport/input/plotting_windows/{point}_plotparams.yaml",
               point = config["points"]),
        expand("synbreaksupport/output/kmer_tables/{point}.pickle",
               point = config["points"]),
        expand("synbreaksupport/output/point_bedfiles/{point}.bed",
               point = config["points"]),
        # get the CCS bamfile
        expand("synbreaksupport/output/bams/{point}_CCS.sorted.bam",
               point = config["points"]),
        # get the CLR bamfile
        expand("synbreaksupport/output/bams/{point}_CLR.sorted.bam",
               point = config["points"]),
        # plot
        expand("synbreaksupport/output/plots/{point}_gapplot.jpg",
               point = config["points"])

rule get_scaffold_sizes:
    input:
        assem    = config["assembly"]
    output:
        scafsize = "synbreaksupport/input/scafsizes.txt"
    threads: 1
    run:
        with open(output.scafsize, "w") as f:
            print("scaf\tscaflength", file = f)
            for record in SeqIO.parse(input.assem, "fasta"):
                print("{}\t{}".format(record.id, len(record.seq)), file = f)

rule get_plotting_info:
    """
    makes a yaml file of all the info we need to plot for this.
    """
    input:
        scafsize = "synbreaksupport/input/scafsizes.txt"
    output:
        plotparams  = "synbreaksupport/input/plotting_windows/{point}_plotparams.yaml"
    threads: 1
    params:
        k = config["kmer"],
    run:
        df = pd.read_csv(input.scafsize, header = 0, sep = "\t")
        print("df")
        print(df)

        scaf1 = wildcards.point.split("-")[0]
        pos1  = int(wildcards.point.split("-")[1])
        scaf2 = wildcards.point.split("-")[2]
        pos2  = int(wildcards.point.split("-")[3])

        scaf1_len = int(df[df["scaf"] == scaf1].squeeze()["scaflength"])
        scaf2_len = int(df[df["scaf"] == scaf2].squeeze()["scaflength"])

        plotdb = {"point":  wildcards.point,
                  "pair" :  "{}-{}".format(scaf1, scaf2),
                  "xaxis":{ "scaf" : scaf1,
                            "min"  : int(max(0, pos1-(config["window"]/2))),
                            "break": pos1,
                            "max"  : int(min(scaf1_len, pos1+(config["window"]/2)))
                          },
                  "yaxis":{ "scaf" : scaf2,
                            "min"  : int(max(0, pos2-(config["window"]/2))),
                            "break": pos2,
                            "max"  : int(min(scaf2_len, pos2+(config["window"]/2)))
                          }
                  }
        print(plotdb)
        plot_params = "synbreaksupport/input/plotting_windows/{}_plotparams.yaml".format(
                     wildcards.point)
        # now save the plotting params since we have all the info we need
        with open(plot_params, "w") as thisfile:
            doc = yaml.dump(plotdb, thisfile)

rule generate_table_from_pairs:
    """
    This just saves a table of kmers that are shared between both pairs
    """
    input:
        assem       = config["assembly"],
        plotparams  = "synbreaksupport/input/plotting_windows/{point}_plotparams.yaml"
    output:
        pointpickle = "synbreaksupport/output/kmer_tables/{point}.pickle"
    threads: 1
    params:
        k = config["kmer"]
    run:
        print("processing {}".format(wildcards.point))
        scaf1 = wildcards.point.split("-")[0]
        pos1  = int(wildcards.point.split("-")[1])
        scaf2 = wildcards.point.split("-")[2]
        pos2  = int(wildcards.point.split("-")[3])

        with open(input.plotparams) as handle:
            plotyaml = yaml.full_load(handle)

        # super slow inplementation, just because it is easier to implement
        seq1 = ""
        seq2 = ""
        for record in SeqIO.parse(input.assem, "fasta"):
            if record.id == scaf1:
                seq1 = str(record.seq)
            elif record.id == scaf2:
                seq2 = str(record.seq)
        kmerdict = {}
        # go through sca1
        print("  - going through scaffold 1")
        for i in range(plotyaml["xaxis"]["min"], plotyaml["xaxis"]["max"]-config["kmer"]+1):
            thiskmer = canonical(seq1[i:i+params.k].upper())
            if thiskmer not in kmerdict:
                kmerdict[thiskmer] = {scaf1: set(), scaf2: set()}
            kmerdict[thiskmer][scaf1].add(i)
        # now go through sca2
        print("  - going through scaffold 2")
        for i in range(plotyaml["yaxis"]["min"], plotyaml["yaxis"]["max"]-config["kmer"]+1):
            thiskmer = canonical(seq2[i:i+params.k].upper())
            if thiskmer in kmerdict:
                kmerdict[thiskmer][scaf2].add(i)
        # now delete the entries that don't have shared kmers
        print("  - finding unnecessary kmers")
        delete_set = set()
        lengthkmerdict = len(kmerdict)
        for key in kmerdict:
            if ( (len(kmerdict[key][scaf1]) == 0) or (len(kmerdict[key][scaf2]) == 0) ):
                delete_set.add(key)
        print("  - deleting unnecessary kmers")
        for delthis in delete_set:
            del kmerdict[delthis]
        # now save this to a file
        print("  - saving the yaml to a file")
        outfile = "synbreaksupport/output/kmer_tables/{}.pickle".format(wildcards.point)

        with open(outfile, "wb") as f:
            pickle.dump(kmerdict, f)

rule generate_bed_from_params:
    input:
        plotparams  = "synbreaksupport/input/plotting_windows/{point}_plotparams.yaml"
    output:
        outbed      = "synbreaksupport/output/point_bedfiles/{point}.bed"
    threads: 1
    run:
        with open(input.plotparams) as handle:
            plotyaml = yaml.full_load(handle)
        # make a bed file
        with open(output.outbed, "w") as handle:
            print("{}\t{}\t{}".format(
                plotyaml["xaxis"]["scaf"],
                plotyaml["xaxis"]["min"],
                plotyaml["xaxis"]["max"]),
                  file = handle)
            print("{}\t{}\t{}".format(
                plotyaml["yaxis"]["scaf"],
                plotyaml["yaxis"]["min"],
                plotyaml["yaxis"]["max"]),
                  file = handle)

rule generate_ccs_bam_from_bed:
    """
    This makes a smaller bam file of the reads.
      Only keeps the best and primary reads.
    """
    input:
        ccsbam      = config["ccs"],
        bed      = "synbreaksupport/output/point_bedfiles/{point}.bed",
        plotparams  = "synbreaksupport/input/plotting_windows/{point}_plotparams.yaml"
    output:
        bam         = "synbreaksupport/output/bams/{point}_CCS.sorted.bam"
    threads: workflow.cores - 1
    shell:
        """
        samtools view -hb -F 2308 -L {input.bed} {input.ccsbam} > {output.bam}
        """

rule index_ccs_bam:
    """
    This makes a smaller bam file of the reads.
      Only keeps the best and primary reads.
    """
    input:
        bam = "synbreaksupport/output/bams/{point}_CCS.sorted.bam"
    output:
        bai = "synbreaksupport/output/bams/{point}_CCS.sorted.bam.bai"
    threads: workflow.cores - 1
    shell:
        """
        samtools index {input.bam}
        """

rule generate_clr_bam_from_bed:
    """
    This makes a smaller bam file of the reads.
      Only keeps the best and primary reads.
    """
    input:
        subbam      = config["subreads"],
        bed      = "synbreaksupport/output/point_bedfiles/{point}.bed",
        plotparams  = "synbreaksupport/input/plotting_windows/{point}_plotparams.yaml"
    output:
        bam         = "synbreaksupport/output/bams/{point}_CLR.sorted.bam"
    threads: workflow.cores - 1
    shell:
        """
        samtools view -hb -F 2308 -L {input.bed} {input.subbam} > {output.bam}
        """

rule index_clr_bam:
    input:
        bam = "synbreaksupport/output/bams/{point}_CLR.sorted.bam"
    output:
        bai = "synbreaksupport/output/bams/{point}_CLR.sorted.bam.bai"
    threads: workflow.cores - 1
    shell:
        """
        samtools index {input.bam}
        """

def plot_reads_in_panel(thisdf, panel):
    import matplotlib.patches as mplpatches
    thisdf = thisdf.reset_index(drop=True)
    done = False
    # plot until we plot em all
    thickness = 0.75
    row_counter = 1
    while not done:
        prev_end = 0
        for index, row in thisdf.iterrows():
            if row["start"] > prev_end:
                plot_bottom = row_counter + (1-thickness)/2
                plot_left = row["start"]
                plot_width = row["end"] - row["start"] + 1
                plot_height = thickness
                plot_color = "#525252" if row["crosses_point"] else "#e2e2e2"
                rectangle1=mplpatches.Rectangle((plot_left,plot_bottom),
                                        plot_width, plot_height,
                                        linewidth=0,
                                        facecolor= plot_color
                                                )
                panel.add_patch(rectangle1)
                thisdf.loc[index, "plotted"] = True
                prev_end = row["end"]
        prev_end = 0
        row_counter += 1
        print("plotted {} reads".format(len(thisdf[thisdf["plotted"] == True])))
        thisdf = thisdf[thisdf["plotted"] == False]
        thisdf = thisdf.reset_index(drop = True)
        if len(thisdf) == 0:
            done = True
    panel.set_ylim([0, 1.15 * row_counter])

def plot_reads_in_panel_y(thisdf, panel):
    import matplotlib.patches as mplpatches
    thisdf = thisdf.reset_index(drop=True)
    done = False
    # plot until we plot em all
    thickness = 0.75
    row_counter = 1
    while not done:
        prev_end = 0
        for index, row in thisdf.iterrows():
            if row["start"] > prev_end:
                plot_left = row_counter + (1-thickness)/2
                plot_bottom = row["start"]
                plot_height = row["end"] - row["start"] + 1
                plot_width = thickness
                plot_color = "#525252" if row["crosses_point"] else "#e2e2e2"
                rectangle1=mplpatches.Rectangle((plot_left,plot_bottom),
                                        plot_width, plot_height,
                                        linewidth=0,
                                        facecolor= plot_color
                                                )
                panel.add_patch(rectangle1)
                thisdf.loc[index, "plotted"] = True
                prev_end = row["end"]
        prev_end = 0
        row_counter += 1
        print("plotted {} reads".format(len(thisdf[thisdf["plotted"] == True])))
        thisdf = thisdf[thisdf["plotted"] == False]
        thisdf = thisdf.reset_index(drop = True)
        if len(thisdf) == 0:
            done = True
    panel.set_xlim([0, 1.15 * row_counter])

rule plot_the_synteny_plot:
    input:
        subscaffold = "synbreaksupport/output/kmer_tables/{point}.pickle",
        plotparams  = "synbreaksupport/input/plotting_windows/{point}_plotparams.yaml",
        ccsbam      = "synbreaksupport/output/bams/{point}_CCS.sorted.bam",
        bai         = "synbreaksupport/output/bams/{point}_CCS.sorted.bam.bai",
        clrbam      = "synbreaksupport/output/bams/{point}_CLR.sorted.bam",
        clrbai      = "synbreaksupport/output/bams/{point}_CLR.sorted.bam.bai",
    output:
        pdf         = "synbreaksupport/output/plots/{point}_gapplot.jpg"
    params:
        k = config["kmer"]
    threads: 1
    run:
        print("plotting {}".format(wildcards.point))
        scaf1 = wildcards.point.split("-")[0]
        pos1  = int(wildcards.point.split("-")[1])
        scaf2 = wildcards.point.split("-")[2]
        pos2  = int(wildcards.point.split("-")[3])

        with open(input.plotparams) as handle:
            plotyaml = yaml.full_load(handle)
        kmerdict = pickle.load( open(input.subscaffold, "rb" ) )

        import pandas as pd
        import seaborn as sns; sns.set()
        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
        from matplotlib.ticker import StrMethodFormatter, NullFormatter
        from matplotlib.patches import Ellipse
        import matplotlib.patches as mplpatches
        import numpy as np
        import pysam
        # set seaborn stuff
        #sns.set(rc={'text.usetex' : True})
        sns.set_style("ticks", {'font.family': ['sans-serif'],
                                    'font.sans-serif': ['Helvetica'],
                                    'grid.color': '.95'})
        # Preserve the vertical order of embedded images:
        matplotlib.rcParams['image.composite_image'] = False
        # text as font in pdf
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

        in_ccs = set()
        # ccsbam
        ccsbam = pysam.AlignmentFile(input.ccsbam, "rb")
        # X AXIS
        plotting_axes = {}
        aln_lens = []
        # get the reads that cross the points
        for thisaxis in ["yaxis", "xaxis"]:
            ccs_that_cross_point = set()
            for read in ccsbam.fetch(
                    plotyaml[thisaxis]["scaf"],
                    int(plotyaml[thisaxis]["break"])-1,
                    int(plotyaml[thisaxis]["break"])+1):
                ccs_that_cross_point.add(read.query_name)
            # now make a dict of x-axis plotting info
            plotting_data = []
            for read in ccsbam.fetch(
                    plotyaml[thisaxis]["scaf"],
                    int(plotyaml[thisaxis]["min"]),
                    int(plotyaml[thisaxis]["max"])):
                in_ccs.add(read.query_name.replace("/ccs", ""))
                plotting_data.append({
                    "readname":  read.query_name,
                    "crosses_point": True if read.query_name in ccs_that_cross_point else False,
                    "positions": find_ranges(read.positions),
                    "strand": "-" if read.is_reverse else "+",
                    "mapq": read.mapq,
                    "start": read.reference_start,
                    "end": read.reference_end,
                    "plotted": False
                    })
                thislen = int(read.reference_end) - int(read.reference_start)
                aln_lens.append(thislen)
            plotting_axes[thisaxis] = pd.DataFrame.from_dict(plotting_data)
            print(plotting_axes[thisaxis])

        mean_ccs_len = np.mean(aln_lens)
        # clrbam
        clrbam = pysam.AlignmentFile(input.clrbam, "rb")
        # X AXIS
        clr_axes = {}
        # get the reads that cross the points
        for thisaxis in ["yaxis", "xaxis"]:
            ccs_that_cross_point = set()
            for read in clrbam.fetch(
                    plotyaml[thisaxis]["scaf"],
                    int(plotyaml[thisaxis]["break"])-1,
                    int(plotyaml[thisaxis]["break"])+1):
                ccs_that_cross_point.add(read.query_name)
            # now make a dict of x-axis plotting info
            plotting_data = []
            for read in clrbam.fetch(
                    plotyaml[thisaxis]["scaf"],
                    int(plotyaml[thisaxis]["min"]),
                    int(plotyaml[thisaxis]["max"])):
                readstem = "/".join(read.query_name.split("/")[0:2])
                if readstem not in in_ccs:
                    thisend   = read.reference_end
                    thisstart = read.reference_start
                    if (thisend - thisstart) >= mean_ccs_len:
                        plotting_data.append({
                            "readname":  read.query_name,
                            "crosses_point": True if read.query_name in ccs_that_cross_point else False,
                            "positions": find_ranges(read.positions),
                            "strand": "-" if read.is_reverse else "+",
                            "mapq": read.mapq,
                            "start": read.reference_start,
                            "end": read.reference_end,
                            "plotted": False
                            })
            clr_axes[thisaxis] = pd.DataFrame.from_dict(plotting_data)
            print(clr_axes[thisaxis])

        # now make a scatter plot
        figWidth = 8
        figHeight = 8
        plt.figure(figsize=(figWidth,figHeight))
        #set the panel dimensions
        panelWidth = 4
        panelHeight = 4
        #find the margins to center the panel in figure
        #leftMargin = (figWidth - panelWidth)/2
        #bottomMargin = ((figHeight - panelHeight)/2)
        leftMargin = 1.5
        bottomMargin = 1.5

        plt.gcf().text((leftMargin + panelWidth + 0.2 + (panelWidth/8) )/figWidth,
             (bottomMargin + panelHeight + 0.2 + (panelHeight/8) )/figHeight,
               "Read\nDepth", horizontalalignment='center', fontsize = 14,
             verticalalignment='center')

        panel1 = plt.axes([leftMargin/figWidth, #left
                             bottomMargin/figHeight,    #bottom
                             panelWidth/figWidth,   #width
                             panelHeight/figHeight])     #height
        panel1.tick_params(axis='both',which='both',
                            bottom=True, labelbottom=True,
                            left=True, labelleft=True,
                            right=True, labelright=False,
                            top=True, labeltop=False)
        panelx = plt.axes([leftMargin/figWidth, #left
                             (bottomMargin + panelHeight + 0.1)/figHeight,    #bottom
                             panelWidth/figWidth,   #width
                             (panelHeight/8)/figHeight])     #height
        panelx.tick_params(axis='both',which='both',
                            bottom=False, labelbottom=False,
                            left=False, labelleft=False,
                            right=True, labelright=True,
                            top=True, labeltop=False)
        panelx2 = plt.axes([leftMargin/figWidth, #left
                             (bottomMargin + panelHeight + 0.1 + (panelHeight/8) + 0.1)/figHeight,    #bottom
                             panelWidth/figWidth,   #width
                             (panelHeight/8)/figHeight])     #height
        panelx2.tick_params(axis='both',which='both',
                            bottom=False, labelbottom=False,
                            left=False, labelleft=False,
                            right=True, labelright=True,
                            top=True, labeltop=True)

        panely = plt.axes([(leftMargin + panelWidth + 0.1)/figWidth, #left
                             bottomMargin/figHeight,    #bottom
                             (panelWidth/8)/figWidth,   #width
                             panelHeight/figHeight])     #height
        panely.tick_params(axis='both',which='both',
                            bottom=False, labelbottom=False,
                            left=False, labelleft=False,
                            right=True, labelright=False,
                            top=True, labeltop=True)
        panely2 = plt.axes([(leftMargin + panelWidth + 0.1 + (panelWidth/8) + 0.1)/figWidth, #left
                             bottomMargin/figHeight,    #bottom
                             (panelWidth/8)/figWidth,   #width
                             panelHeight/figHeight])     #height
        panely2.tick_params(axis='both',which='both',
                            bottom=False, labelbottom=False,
                            left=False, labelleft=False,
                            right=True, labelright=True,
                            top=True, labeltop=True)

        panel1.set_xlim([plotyaml["xaxis"]["min"], plotyaml["xaxis"]["max"]])
        panel1.set_ylim([plotyaml["yaxis"]["min"], plotyaml["yaxis"]["max"]])
        panel1.invert_yaxis()

        numElems = 9
        xidx = np.linspace(plotyaml["xaxis"]["min"], plotyaml["xaxis"]["max"], numElems)
        yidx = np.linspace(plotyaml["yaxis"]["min"], plotyaml["yaxis"]["max"], numElems)
        xlabels = [round(x/1000000, 3) for x in xidx]
        ylabels = [round(y/1000000, 3) for y in yidx]

        #panel1.xaxis.set_label_position("top")
        panel1.set_xticks(xidx)
        panel1.set_xticklabels(xlabels, fontsize=8)
        panel1.set_xlabel("{} Mb".format(plotyaml["xaxis"]["scaf"]))

        panel1.set_yticks(yidx)
        panel1.set_yticklabels(ylabels, fontsize=8)
        panel1.set_ylabel("{} Mb".format(plotyaml["yaxis"]["scaf"]))

        # plot the break
        panel1.axvline(x=plotyaml["xaxis"]["break"], color = "#b53630", ls = "--", alpha = 0.5)
        panel1.axhline(y=plotyaml["yaxis"]["break"], color = "#b53630", ls = "--", alpha = 0.5)
        xthou = (plotyaml["xaxis"]["max"] - plotyaml["xaxis"]["min"])/500
        ythou = (plotyaml["yaxis"]["max"] - plotyaml["yaxis"]["min"])/500

        xscaf = plotyaml["xaxis"]["scaf"]
        yscaf = plotyaml["yaxis"]["scaf"]
        counter = 1
        lendict = len(kmerdict)
        for thiskey in kmerdict:
            if len(kmerdict[thiskey][xscaf]) <= 100:
                for xpos in kmerdict[thiskey][xscaf]:
                   if len(kmerdict[thiskey][yscaf]) <= 100:
                       for ypos in kmerdict[thiskey][yscaf]:
                           thisEllipse = Ellipse(xy = [xpos, ypos],
                                                 width = xthou,
                                                 height = ythou,
                                                 angle = 0,
                                                 alpha = 0.25,
                                                 linewidth = 0
                                                )
                           panel1.add_artist(thisEllipse)
            print("    done with kmer {} of {}                      ".format(counter, lendict), end = "\r")
            counter += 1
        print("    done with kmer {} of {}                      ".format(counter, lendict), end = "\r")
        print()

        # now plot the x-axis reads
        print("plotting CCS x-axis reads")
        panelx.set_xlim([plotyaml["xaxis"]["min"], plotyaml["xaxis"]["max"]])
        panelx.tick_params(axis = "y", which = 'major', labelsize = 6)
        panelx.set_xticks(xidx)
        plot_reads_in_panel(plotting_axes["xaxis"].sort_values("start"), panelx)
        panelx.yaxis.set_label_position("left")
        panelx.set_ylabel("CCS\nReads", rotation=90, fontsize=8)

        # now plot the x-axis clr reads
        print("plotting CLR x-axis reads")
        panelx2.set_xlim([plotyaml["xaxis"]["min"], plotyaml["xaxis"]["max"]])
        panelx2.tick_params(axis = "y", which = 'major', labelsize = 6)
        panelx2.xaxis.set_label_position("top")
        panelx2.set_xticks(xidx)
        panelx2.set_xticklabels(xlabels, fontsize=8)
        panelx2.set_xlabel("{} Mb".format(plotyaml["xaxis"]["scaf"]))
        plot_reads_in_panel(clr_axes["xaxis"].sort_values("start"), panelx2)
        panelx2.yaxis.set_label_position("left")
        panelx2.set_ylabel("CLR\nReads", rotation=90, fontsize=8)

        # now plot the y-axis reads
        print("Plotting the y-axis reads")
        panely.set_ylim([plotyaml["yaxis"]["min"], plotyaml["yaxis"]["max"]])
        panely.yaxis.set_label_position("right")
        panely.tick_params(axis = "x", which = 'major', labelsize = 6)
        panely.set_yticks(yidx)
        plot_reads_in_panel_y(plotting_axes["yaxis"].sort_values("start"), panely)
        panely.invert_yaxis()
        panely.xaxis.set_label_position("bottom")
        panely.set_xlabel("CCS\nReads", rotation=0, fontsize=8)

        # now plot the y-axis CLR
        print("plotting CLR y-axis reads")
        panely2.set_ylim([plotyaml["yaxis"]["min"], plotyaml["yaxis"]["max"]])
        panely2.yaxis.set_label_position("right")
        panely2.tick_params(axis = "x", which = 'major', labelsize = 6)
        panely2.set_yticks(yidx)
        panely2.set_yticklabels(ylabels, fontsize=8)
        panely2.set_ylabel("{} Mb".format(plotyaml["yaxis"]["scaf"]))
        plot_reads_in_panel_y(clr_axes["yaxis"].sort_values("start"), panely2)
        panely2.invert_yaxis()
        panely2.xaxis.set_label_position("bottom")
        panely2.set_xlabel("CLR\nReads", rotation=0, fontsize=8)

        print("saving {}".format(output.pdf))
        plt.savefig(output.pdf, dpi = 600)
        print("saved {}".format(output.pdf))

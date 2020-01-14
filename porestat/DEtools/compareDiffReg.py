
import argparse
import sys
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from collections import Counter
import numpy as np
from matplotlib import ticker
import pandas as pd
import matplotlib.pyplot as plt

def printParallelCoords(df, cols, outname, plotWidth, plotLength, invRank=True, fcRange=None):

    x = [i for i, _ in enumerate(cols)]
    # create dict of categories: colours

    # Create (X-1) sublots along x axis
    fig, axes = plt.subplots(1, len(x) - 1, sharey=True, figsize=(plotWidth, plotLength))

    # Get min, max and range for each column
    # Normalize the data for each column
    min_max_range = {}

    for col in cols:

        if "log2FC" in col and fcRange != None:
            colMin = fcRange[0]
            colMax = fcRange[1]
            colRange = colMax-colMin

            min_max_range[col] = [colMin, colMax, colRange]
            df[col] = np.true_divide(df[col] - colMin, colRange)
        else:

            colMin = df[col].min()
            colMax = df[col].max()

            min_max_range[col] = [colMin, colMax, np.ptp(df[col])]
            df[col] = np.true_divide(df[col] - colMin, np.ptp(df[col]))

        print(col, df[col].min(), df[col].max())


    # Plot each row
    for i, ax in enumerate(axes):
        for idx in df.index:
            ax.plot(x, df.loc[idx, cols], c='#1f77b4')
        ax.set_xlim([x[i], x[i + 1]])


    # Set the tick positions and labels on y axis for each plot
    # Tick positions based on normalised data
    # Tick labels are based on original data
    def set_ticks_for_axis(dim, ax, ticks, useFCRange):
        min_val, max_val, val_range = min_max_range[cols[dim]]

        if "log2FC" in cols[dim]:
            min_val = fcRange[0]
            max_val = fcRange[1]
            val_range = fcRange[1]-fcRange[0]

        step = val_range / float(ticks - 1)
        tick_labels = [round(min_val + step * i, 2) for i in range(ticks)]

        norm_min = df[cols[dim]].min()
        norm_range = np.ptp(df[cols[dim]])

        if "log2FC" in cols[dim]:
            norm_min = 0.0
            norm_range = 1.0

        print(cols[dim], min_val, max_val, val_range, norm_min, norm_range)

        norm_step = norm_range / float(ticks - 1)
        ticks = [round(norm_min + norm_step * i, 2) for i in range(ticks)]
        ax.yaxis.set_ticks(ticks)
        ax.set_yticklabels(tick_labels)

        print(cols[dim], ticks)
        print(cols[dim], tick_labels)


    for dim, ax in enumerate(axes):
        ax.xaxis.set_major_locator(ticker.FixedLocator([dim]))

        ufc = fcRange != None and "log2FC" in cols[dim]

        set_ticks_for_axis(dim, ax, ticks=10, useFCRange=ufc)
        ax.set_xticklabels([cols[dim]])

    # Move the final axis' ticks to the right-hand side
    ax = plt.twinx(axes[-1])
    dim = len(axes)
    ax.xaxis.set_major_locator(ticker.FixedLocator([x[-2], x[-1]]))
    set_ticks_for_axis(dim, ax, ticks=10, useFCRange=fcRange != None and "log2FC" in cols[-1])
    ax.set_xticklabels([cols[-2], cols[-1]])

    #for i, ax in enumerate(axes):
    #    if "rank" in cols[i]:
    #        ax.invert_yaxis()

    # Remove space between subplots
    plt.subplots_adjust(wspace=0)

    # Add legend to plot
    #plt.legend()

    plt.grid(b=None)
    plt.savefig(outname, bbox_inches="tight")
    plt.close()

    print(outname)

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d1', '--de1', type=argparse.FileType('r'), required=True, help='de files')
    parser.add_argument('-d2', '--de2', type=argparse.FileType('r'), required=True, help='de files')
    parser.add_argument('-o', '--output', type=str, required=True, help='de files')

    parser.add_argument("-m", "--method", type=str, required=True)
    parser.add_argument("-n1", "--df1name", type=str, required=False, default="DF1")
    parser.add_argument("-n2", "--df2name", type=str, required=False, default="DF2")

    parser.add_argument("-c", "--common", required=False, default=False, action='store_true')
    parser.add_argument("-n", "--num", type=int, required=False, default=50)

    args = parser.parse_args()

    indf1 = DataFrame.parseFromFile(args.de1.name, skipChar='#', replacements = {
        "None": None,
        "": None,
        "NA": None
    })

    indf2 = DataFrame.parseFromFile(args.de2.name, skipChar='#', replacements = {
        "None": None,
        "": None,
        "NA": None
    })



    df1Cols = indf1.getHeader()
    df2Cols = indf2.getHeader()

    df1SpecialCols = [x for x in df1Cols if any(["PVAL" in x, "FC" in x])]
    df2SpecialCols = [x for x in df2Cols if any(["PVAL" in x, "FC" in x])]


    for x in df1SpecialCols:
        xn = x.split("_")
        xn.insert(1, "DF1")
        xn = "_".join(xn)


    for x in df2SpecialCols:
        xn = x.split("_")
        xn.insert(1, "DF2")
        xn = "_".join(xn)

    geneIdCol = "id"

    def loadDF(indf):
        gene2fc = {}
        allgenes = []
        for row in indf:
            geneName = row[geneIdCol]
            geneFC = float(row[args.method + "_log2FC"])
            genePVal = float(row[args.method + "_ADJ.PVAL"])

            allgenes.append((geneName, geneFC, genePVal))
            gene2fc[geneName] = geneFC

        return allgenes, gene2fc


    allGenesDF1, gene2fcDF1 = loadDF(indf1)
    allGenesDF2, gene2fcDF2 = loadDF(indf2)


    def getGenes(allGenes):

        #filter significants
        allGenes = [x for x in allGenes if x[2] < 0.05]

        #filter down and up
        upGenes = [x for x in allGenes if x[1] >= 0]
        downGenes = [x for x in allGenes if x[1] <= 0]

        return upGenes, downGenes

    selUpDF1, selDownDF1 = getGenes(allGenesDF1)
    selUpDF2, selDownDF2 = getGenes(allGenesDF2)

    selUpDF1 = sorted(selUpDF1, key=lambda x: (x[2], x[1]))
    selUpDF2 = sorted(selUpDF2, key=lambda x: (x[2], x[1]))

    selDownDF1 = sorted(selDownDF1, key=lambda x: (x[2], x[1]))
    selDownDF2 = sorted(selDownDF2, key=lambda x: (x[2], x[1]))



    def selElems(inlist, num):

        if len(inlist) < num or num == -1:
            return [x for x in inlist]

        return inlist[0:num]


    def makeIntersectStat(selDF1, selDF2):

        stat = []

        useElems = min([len(selDF1), len(selDF2)])
        for i in range(1, useElems+1):

            df1Sel = selElems(selDF1, i)
            df2Sel = selElems(selDF2, i)

            df1SelNames = [x[0] for x in df1Sel]
            df2SelNames = [x[0] for x in df2Sel]

            intersectCount = len(set(df1SelNames).intersection(df2SelNames))
            unionCount = len(set(df1SelNames).union(df2SelNames))

            stat.append( (i, intersectCount/unionCount, intersectCount/i ) )

        return stat

    statUP = makeIntersectStat(selUpDF1, selUpDF2)
    statDOWN = makeIntersectStat(selDownDF1, selDownDF2)


    plt.figure()
    plt.plot([x[0] for x in statUP], [x[1] for x in statUP], label="UP union")
    plt.plot([x[0] for x in statUP], [x[2] for x in statUP], label="UP single")

    plt.plot([x[0] for x in statDOWN], [x[1] for x in statDOWN], label="DOWN union")
    plt.plot([x[0] for x in statDOWN], [x[2] for x in statDOWN], label="DOWN single")
    plt.legend()
    plt.show()
    plt.close()

    exit()




    numSelected = args.num

    selUpDF1 = selElems(selUpDF1, numSelected)
    selUpDF2 = selElems(selUpDF2, numSelected)

    selDownDF1 = selElems(selDownDF1, numSelected)
    selDownDF2 = selElems(selDownDF2, numSelected)

    allGeneNamesUp = None
    allGeneNamesDown = None

    if args.common:
        allGeneNamesUp = set([x[0] for x in selUpDF1]).intersection([x[0] for x in selUpDF2])
        allGeneNamesDown = set([x[0] for x in selDownDF1]).intersection([x[0] for x in selDownDF2])
    else:
        allGeneNamesUp = set([x[0] for x in selUpDF1] + [x[0] for x in selUpDF2])
        allGeneNamesDown = set([x[0] for x in selDownDF1] + [x[0] for x in selDownDF2])

    geneUp2Ranks = []
    geneDown2Ranks = []

    selUpDF1Names = [x[0] for x in selUpDF1]
    selUpDF2Names = [x[0] for x in selUpDF2]
    selDownDF1Names = [x[0] for x in selDownDF1]
    selDownDF2Names = [x[0] for x in selDownDF2]

    for x in allGeneNamesUp:
        df1idx = numSelected+1
        df2idx = numSelected+1

        try:
            df1idx = selUpDF1Names.index(x)
        except:
            pass

        try:
            df2idx = selUpDF2Names.index(x)
        except:
            pass

        geneUp2Ranks.append(
            (x, df1idx, df2idx )
        )

    for x in allGeneNamesDown:
        df1idx = numSelected+1
        df2idx = numSelected+1

        try:
            df1idx = selDownDF1Names.index(x)
        except:
            pass

        try:
            df2idx = selDownDF2Names.index(x)
        except:
            pass

        geneDown2Ranks.append(
            (x, df1idx, df2idx )
        )

    geneUp2Ranks = sorted(geneUp2Ranks, key=lambda x: (x[1], x[2]))
    geneDown2Ranks = sorted(geneDown2Ranks, key=lambda x: (x[1], x[2]))

    #for x in geneUp2Ranks:
    #    print(x)

    #print()
    #for x in geneDown2Ranks:
    #    print(x)

    updata = []


    dfnames = [args.df1name + " log2FC", args.df1name +' rank', args.df2name +' rank', args.df2name + " log2FC"]

    for x in geneUp2Ranks:
        updata.append({
            args.df1name + " log2FC": gene2fcDF1[x[0]],
            args.df1name +' rank': x[1],
            args.df2name +' rank': x[2],
            args.df2name + " log2FC": gene2fcDF2[x[0]],
            'name': "Top DE Up Genes (n={} of {})".format(args.num, len(allGeneNamesUp)),
        })

    downdata = []
    for x in geneDown2Ranks:
        downdata.append({
            args.df1name + " log2FC": gene2fcDF1[x[0]],
            args.df1name +' rank': x[1],
            args.df2name +' rank': x[2],
            args.df2name + " log2FC": gene2fcDF2[x[0]],
            'name': "Top DE Down Genes (n={} of {})".format(args.num, len(allGeneNamesDown)),
        })



    plotLength = 12

    if args.num == -1:
        plotLength = 36


    df = pd.DataFrame(updata)
    df.to_csv(args.output + ".up.tsv")

    df = pd.DataFrame(downdata)
    df.to_csv(args.output + ".down.tsv")

    if args.num == -1:

        upRankDiffs = []
        for x in geneUp2Ranks:
            rankDiff = abs(x[1]-x[2])
            upRankDiffs.append(rankDiff)

        downRankDiffs = []
        for x in geneDown2Ranks:
            rankDiff = abs(x[1] - x[2])
            downRankDiffs.append(rankDiff)

        fig = plt.figure()
        plt.hist(upRankDiffs, bins=len(upRankDiffs), label="UP rank diffs (n={})".format(len(upRankDiffs)), normed=True, histtype="step", cumulative=True)
        plt.hist(downRankDiffs, bins=len(downRankDiffs), label="DOWN rank diffs (n={})".format(len(upRankDiffs)), normed=True, histtype="step", cumulative=True)
        plt.legend()
        plt.savefig(args.output + ".rankdiff.png", bbox_inches="tight")
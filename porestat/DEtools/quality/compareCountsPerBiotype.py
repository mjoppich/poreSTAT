
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import matplotlib as mpl
import argparse
import math
import natsort

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

# mpl.style.use("seaborn")


def makeConditionNames(conditions, inHeaders):

    allSamples = []
    for x in conditions:
        for y in x:
            if not y in allSamples:
                allSamples.append(y)


    newconds = args.conditions

    dfCols = indf.getHeader()

    allDerivedSamples = []
    replaceSamples = {}
    for sample in allSamples:

        if sample in inHeaders:
            replaceSamples[sample] = sample
            allDerivedSamples.append(sample)

        else:
            sampleAdded = False
            for colName in inHeaders:
                if colName.endswith(sample[1:]):
                    sampleAdded = True
                    replaceSamples[sample] = colName
                    allDerivedSamples.append(colName)
                    break
                elif colName.replace("//", "/").endswith(sample.replace("//", "/")[1:]):
                    sampleAdded = True
                    replaceSamples[sample] = colName
                    allDerivedSamples.append(colName)
                    break

            if not sampleAdded:
                print("Could not assign sample", sample)
                print(inHeaders)

                if not args.ignoreMissing:
                    exit(-1)

    allSamples = allDerivedSamples

    for sample in replaceSamples:
        print(sample, "-->", replaceSamples[sample])

    dfGroups = args.conditions
    newgroups = []
    for gelem in dfGroups:
        if gelem in replaceSamples:  # may only not be included, if not added -> error before
            newgroups.append(replaceSamples[gelem])

    return newgroups


def autolabel(ax, rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()

        ltext = str(round(height, 2))

        if height < 0.01:
            ltext = ""

        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                ltext,
                ha='left', va='bottom', rotation=45)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-d', '--conditions' ,nargs='+', type=str,  required=True, help='alignment files')
    parser.add_argument('-l', '--last', action="store_true", default=False)
    parser.add_argument('-i', '--ignoreMissing', action="store_true", default=False)
    parser.add_argument('-o', '--output', nargs='+', type=str, required=False, help="output base")
    parser.add_argument("-b", '--biotype', type=argparse.FileType('r'), required=True)

    parser.add_argument('-p', '--pdf', action="store_true", default=False)

    args = parser.parse_args()

    if args.output == None:
        args.output = [counts.name for counts in args.counts]

    geneID2Type = {}
    allTypes = set()

    for lidx, line in enumerate(args.biotype):

        if lidx == 0:
            continue

        aline = line.strip().split("\t")

        ensID = aline[0]
        geneSym = aline[1]
        geneType = aline[2].strip()

        allTypes.add(geneType)

        if not ensID in geneID2Type:
            geneID2Type[ensID] = geneType

        geneID2Type[geneSym] = geneType

    allTypes = sorted(allTypes)

    for fidx, defile in enumerate(args.counts):
        indf = DataFrame.parseFromFile(defile.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })


        inHeaders = indf.getHeader()
        newconds = makeConditionNames([args.conditions], inHeaders)

        
        cond2allcounts = {}
        cond2ProtCodingCounts = {}


        for cidx, condition in enumerate(sorted(newconds)):
            geneTypeCounts = Counter()

            for row in indf:

                geneCount = row[condition]

                geneSym = None

                if "gene_symbol" in inHeaders:
                    geneSym = row["gene_symbol"]
                elif "Geneid" in inHeaders:
                    geneSym = row["Geneid"]
                elif "id" in inHeaders:
                    geneSym = row["id"]


                if geneSym == None:
                    continue

                geneSymType = geneID2Type.get(geneSym, "Unknown")

                geneTypeCounts[geneSymType] += geneCount

            plotFigSize=(14, 6)

            fig, ax = plt.subplots(figsize=plotFigSize, dpi=80)

            allValues = [geneTypeCounts.get(gtype, 0) for gtype in allTypes]
            allValueCount = sum(allValues)
            ind = range(0, len(allTypes))

            # save for later
            cond2allcounts[condition] = allValueCount
            cond2ProtCodingCounts[condition] = geneTypeCounts.get("protein_coding", 0)


            rects = ax.bar(ind, allValues)
            autolabel(ax, rects)

            plt.title(condition)
            plt.xticks(ind, allTypes, rotation=90)

            plt.xlabel("Gene Type")
            plt.ylabel("Counters of all Genes")


            plt.savefig(args.output[fidx] + ".cpergenes." +str(fidx ) + "." + str(cidx) + ".png", bbox_inches ="tight")
            plt.close()



            fig, ax = plt.subplots(figsize=plotFigSize, dpi=80)

            allValues = [geneTypeCounts.get(gtype, 0)/allValueCount for gtype in allTypes]
            ind = range(0, len(allTypes))

            rects = ax.bar(ind, allValues)
            autolabel(ax, rects)

            plt.title(condition)
            plt.xticks(ind, allTypes, rotation=90)

            plt.xlabel("Gene Type")
            plt.ylabel("Counters of all Genes")


            plt.savefig(args.output[fidx] + ".cpergenes.rel." +str(fidx ) + "." + str(cidx) + ".png", bbox_inches ="tight")
            plt.close()



        fig, ax = plt.subplots(figsize=plotFigSize, dpi=80)

        sortedConds = natsort.natsorted(cond2allcounts)
        allValues = [cond2allcounts[cond] for cond in sortedConds]
        ind = range(0, len(sortedConds))

        rects = ax.bar(ind, allValues)
        autolabel(ax, rects)

        plt.title("Library Size/All Read Counts per Condition")
        plt.xticks(ind, sortedConds, rotation=90)

        plt.xlabel("Condition")
        plt.ylabel("Library Size/All Read Counts")


        plt.savefig(args.output[fidx] + ".cpergenes.all." +str(fidx ) + ".png", bbox_inches ="tight")
        if args.pdf:
            plt.savefig(args.output[fidx] + ".cpergenes.all." +str(fidx ) + ".pdf", bbox_inches ="tight")
        plt.close()



        fig, ax = plt.subplots(figsize=plotFigSize, dpi=80)
        sortedConds = natsort.natsorted(cond2ProtCodingCounts)
        allValues = [cond2ProtCodingCounts[cond] for cond in sortedConds]
        ind = range(0, len(sortedConds))

        rects = ax.bar(ind, allValues)
        autolabel(ax, rects)

        plt.title("Library Size/All Read Counts (protein_coding) per Condition")
        plt.xticks(ind, sortedConds, rotation=90)

        plt.xlabel("Condition")
        plt.ylabel("Library Size/All Read Counts (protein_coding)")


        plt.savefig(args.output[fidx] + ".cpergenes.all_prot." +str(fidx ) + ".png", bbox_inches ="tight")

        if args.pdf:
            plt.savefig(args.output[fidx] + ".cpergenes.all_prot." +str(fidx ) + ".pdf", bbox_inches ="tight")
        plt.close()
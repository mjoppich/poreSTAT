
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import matplotlib as mpl
import argparse
import math

import sys, os
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

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
    for group in dfGroups:
        ngroup = []
        for gelem in group:
            if gelem in replaceSamples:  # may only not be included, if not added -> error before
                ngroup.append(replaceSamples[gelem])

        newgroups.append(ngroup)

    newconds = newgroups

    return newconds


def autolabel(ax, rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                str(round(height, 2)),
                ha='center', va='bottom')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-d', '--conditions', nargs='+', type=str, action='append',  required=True, help='alignment files')
    parser.add_argument('-l', '--last', action="store_true", default=False)
    parser.add_argument('-i', '--ignoreMissing', action="store_true", default=False)
    parser.add_argument('-o', '--output', nargs='+', type=str, required=False, help="output base")
    parser.add_argument("-b", '--biotype', type=argparse.FileType('r'), required=True)

    args = parser.parse_args()

    if args.output == None:
        args.output = [counts.name for counts in args.counts]

    geneID2Type = {}
    allTypes = set()

    for line in args.biotype:

        aline = line.strip().split("\t")

        ensID = aline[0]
        geneSym = aline[1]
        geneType = aline[2]

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
        newconds = makeConditionNames(args.conditions, inHeaders)

        for condition in newconds:
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



            fig, ax = plt.subplots()

            allValues = [geneTypeCounts.get(gtype, 0) for gtype in geneTypeCounts]
            ind = range(0, len(allTypes))

            rects = ax.bar(ind, allValues)
            autolabel(ax, rects)

            plt.title(condition)

            plt.xticks(ind, allTypes, rotation=70)

            plt.xlabel("Category Rel Count")
            plt.ylabel("Assigned Rel Count")


            plt.savefig(args.output[fidx] + ".cpergenes. " +str(fidx ) + ".png", bbox_inches ="tight")
            plt.close()




import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import matplotlib as mpl
import argparse
import math

import sys, os

from upsetplot import from_contents, plot

sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
import venn

if __name__ == '__main__':

    allowedStats = ["ROB_ADJ.PVAL","ROB_log2FC","ROB_log2FC_SIG"]

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d', '--detable', nargs='+', type=argparse.FileType('r'), required=True, help='DE files')
    parser.add_argument('-name', '--name', nargs='+', type=str, required=True, help="name for DE files")
    parser.add_argument('-n', '--top_n', type=int, nargs='+', required=False, default=100)
    parser.add_argument('-o', '--output', type=str, required=False, help="output base")
    parser.add_argument('-s', '--stats', type=str, nargs="+", required=False, default=allowedStats)

    args = parser.parse_args()

    if len(args.detable) != len(args.name):
        print("Mismatch in detable and name lenghts")
        exit(-1)

    for x in args.stats:
        if not x in allowedStats:
            print("Allowed stat values", allowedStats)
            exit(-1)


    foundRes = defaultdict(lambda: defaultdict(list))

    for didx, detableFile in enumerate(args.detable):

        indf = DataFrame.parseFromFile(detableFile.name, skipChar='#', replacements={
            "None": None,
            "": None,
            "NA": None
        })

        print(detableFile.name)

        # try to find ID column
        idColumn = None
        inHeaders = indf.getHeader()

        if "id" in inHeaders:
            idColumn = "id"
        elif "gene_symbol" in inHeaders:
            idColumn = "gene_symbol"


        if idColumn == None:
            continue

        compMethod2Column = {
            'ROB_log2FC': 'ROB_log2FC',
            'ROB_ADJ.PVAL': 'ROB_ADJ.PVAL',
            'ROB_log2FC_SIG': 'ROB_log2FC'
        }


        for topN in args.top_n:
            for compMethod in args.stats:

                print(detableFile.name, compMethod, topN)

                valColumn = compMethod2Column[compMethod]

                topNIDs = []
                for ridx, row in enumerate(indf):

                    eid = row[idColumn]

                    try:
                        eValue = float(row[valColumn])

                        if valColumn in ["ROB_log2FC"]:
                            eValue = abs(eValue)

                        ePValue = 1
                        if compMethod in ["ROB_log2FC_SIG"]:
                            ePValue = float(row["ROB_ADJ.PVAL"])

                            if ePValue > 0.05:
                                continue

                        if row["id"] == "ENSG00000120885":
                            print(detableFile.name, row["id"], eValue, ePValue)

                        topNIDs.append((eid, eValue))

                    except:
                        continue


                sortReversed = valColumn in ["ROB_log2FC"]

                topNIDs = sorted(topNIDs, key=lambda x: x[1], reverse=sortReversed)
                lqCount = len(topNIDs)
                if topN > -1 and len(topNIDs) >= topN:

                    topNElement = topNIDs[topN-1]
                    print("Top N Element Score", topNElement)

                    if sortReversed:
                        lqCount = sum([1 for x in topNIDs if x[1] >= topNElement[1]])
                    else:
                        lqCount = sum([1 for x in topNIDs if x[1] <= topNElement[1]])

                    topNIDs = topNIDs[:topN]

                print(detableFile.name, args.name[didx], compMethod, valColumn, topN, lqCount)

                topNID = topN

                if topN == -1:
                    topNID = "all"

                foundRes[topNID][compMethod].append((args.name[didx], topNIDs, lqCount))

    print("Starting Plotting")
    for topN in foundRes:
        print("TopN", topN, "methods", [x for x in foundRes[topN]])
        for compMethod in foundRes[topN]:

            print(topN, compMethod, "before sets")

            inputSets = [set([y[0] for y in x[1]]) for x in foundRes[topN][compMethod]]

            print(topN, compMethod, "after sets")

            method2genes = {}
            for x in foundRes[topN][compMethod]:
                method2genes[x[0]] = set([y[0] for y in x[1]])
                print(x[0], len(method2genes[x[0]]))

            #print(set(method2genes["RobustDE+Robust"]).difference(method2genes["combined+msEmpiRe_DESeq2"]))

            upIn = from_contents(method2genes)

            print(topN, compMethod, "after content")

            lvls = set([x for x in upIn.index.levshape])
            print("lvls", lvls)

            if len(lvls) <= 2:

                plt.figure()
                plt.title("no data to plot - maybe only 1 or 2 groups?")
                outname = args.output + "." + str(topN) + "." + compMethod

                if not outname.endswith(".png"):
                    outname += ".png"

                plt.savefig(outname, bbox_inches="tight")
                print("finished plot", outname)
                plt.close()
                continue


            plot(upIn,show_counts=True, subset_size="auto")

            print(topN, compMethod, "after plot")

            plt.title("Method: {} (Top {})".format(compMethod, topN))

            outname = args.output + "." + str(topN) + "." + compMethod

            if not outname.endswith(".png"):
                outname += ".png"

            plt.savefig(outname, bbox_inches="tight")
            print("finished plot", outname)
            plt.close()


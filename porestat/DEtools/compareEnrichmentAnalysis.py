

import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import matplotlib as mpl
import argparse
import math

import sys, os
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
import venn

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-p', '--pathways', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-n', '--top_n', nargs='+', type=int, required=False, default=100)
    parser.add_argument('-o', '--output', type=str, required=False, help="output base")


    args = parser.parse_args()


    if len(args.pathways) <= 2 or len(args.pathways) > 5:
        exit()



    pwcount2function = {
        2: venn.venn2,
        3: venn.venn3,
        4: venn.venn4,
        5: venn.venn5,
        6: venn.venn6,
    }


    foundRes = defaultdict(lambda: list())

    for pathwaysFile in args.pathways:
        indf = DataFrame.parseFromFile(pathwaysFile.name, skipChar='#', replacements={
            "None": None,
            "": None,
            "NA": None
        })



        # try to find ID column

        idColumn = None
        inHeaders = indf.getHeader()

        for header in inHeaders:

            if "ID" in header:
                idColumn = header
                break


        if idColumn == None:
            continue

        topNIDs = []
        for ridx, row in enumerate(indf):

            qvalue = float(row.get("qvalue", row.get("qvalues")))

            topNIDs.append((row[idColumn], qvalue))

        topNIDs = sorted(topNIDs, key=lambda x: x[1])

        for topN in args.top_n:

            localTopNIDs = [x for x in topNIDs]

            lqCount = len(localTopNIDs)
            if len(localTopNIDs) >= topN:

                topNElement = localTopNIDs[topN - 1]

                lqCount = sum([1 for x in localTopNIDs if x[1] <= topNElement[1]])

                localTopNIDs = localTopNIDs[:topN]

            print(pathwaysFile.name, topN, len(localTopNIDs))
            foundRes[topN].append((os.path.basename(pathwaysFile.name), localTopNIDs, lqCount))




    for topN in args.top_n:
        vennLabels = venn.generate_petal_labels([set([y[0] for y in x[1]]) for x in foundRes[topN]])
        fig, ax = pwcount2function[len(foundRes[topN])](vennLabels, names=["{fn} (lq={lqc})".format(fn=x[0], lqc=x[2]) for x in foundRes[topN]])

        plt.suptitle("Overlaps for topN={} pathways (by qvalue)".format(topN))

        outname = args.output + "." + str(topN)

        if not outname.endswith(".png"):
            outname += ".png"

        print(outname)
        plt.savefig(outname, bbox_inches ="tight")


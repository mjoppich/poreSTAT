

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
    parser.add_argument('-d', '--detable', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-n', '--top_n', type=int, nargs='+', required=False, default=100)
    parser.add_argument('-o', '--output', type=str, required=False, help="output base")


    args = parser.parse_args()


    if len(args.detable) <= 2 or len(args.detable) > 5:
        exit()



    pwcount2function = {
        2: venn.venn2,
        3: venn.venn3,
        4: venn.venn4,
        5: venn.venn5,
        6: venn.venn6,
    }


    foundRes = defaultdict(lambda: defaultdict(list))

    for detableFile in args.detable:

        if "_raw" in detableFile.name:
            continue


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
            for compMethod in ["ROB_log2FC", "ROB_ADJ.PVAL", "ROB_log2FC_SIG"]:

                print(detableFile.name, compMethod, topN)

                valColumn = compMethod2Column[compMethod]

                topNIDs = []
                for ridx, row in enumerate(indf):

                    eid = row[idColumn]

                    try:
                        eValue = float(row[valColumn])

                        if valColumn in ["ROB_log2FC"]:
                            eValue = abs(eValue)

                        if compMethod in ["ROB_log2FC_SIG"]:
                            ePValue = float(row["ROB_ADJ.PVAL"])

                            if ePValue > 0.05:
                                continue

                        topNIDs.append((eid, eValue))

                    except:
                        continue


                sortReversed = valColumn in ["ROB_log2FC"]

                topNIDs = sorted(topNIDs, key=lambda x: x[1], reverse=sortReversed)

                if len(topNIDs) >= topN:

                    topNElement = topNIDs[topN-1]

                    if sortReversed:
                        lqCount = sum([1 for x in topNIDs if x[1] >= topNElement[1]])
                    else:
                        lqCount = sum([1 for x in topNIDs if x[1] <= topNElement[1]])

                    topNIDs = topNIDs[:topN]

                print(detableFile.name, compMethod, valColumn, topN, lqCount)


                foundRes[topN][compMethod].append((detableFile.name, topNIDs, lqCount))

    for topN in args.top_n:
        for compMethod in foundRes[topN]:

            inputSets = [set([y[0] for y in x[1]]) for x in foundRes[topN][compMethod]]

            print(topN, compMethod)
            #print(set.intersection(*inputSets))

            vennLabels = venn.generate_petal_labels(inputSets)
            fig, ax = pwcount2function[len(foundRes[topN][compMethod])](vennLabels, names=["{fn} (lq={lqc})".format(fn=x[0], lqc=x[2]) for x in foundRes[topN][compMethod]])
            plt.suptitle("Overlaps for topN={} DE genes (by {})".format(topN, compMethod))

            outname = args.output + "." + str(topN) + "." + compMethod

            if not outname.endswith(".png"):
                outname += ".png"

            plt.savefig(outname, bbox_inches ="tight")


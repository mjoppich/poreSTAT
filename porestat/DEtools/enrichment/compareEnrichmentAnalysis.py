

import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import matplotlib as mpl
import argparse
import math

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")


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
        1: venn.venn2,
        2: venn.venn2,
        3: venn.venn3,
        4: venn.venn4,
        5: venn.venn5,
        6: venn.venn6,
    }


    foundRes = defaultdict(lambda: list())
    input2allres = {}

    for pathwaysFile in args.pathways:
        print(pathwaysFile.name)

        indf = DataFrame.parseFromFile(pathwaysFile.name, skipChar='#', replacements={
            "None": None,
            "": None,
            "NA": None
        })

        if indf == None:
            continue



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

            qvaluestr = row.get("qvalue", row.get("qvalues", None))
            qvalue = 1

            if qvaluestr != None and qvaluestr != "NA":
                qvalue = float(qvaluestr)

            topNIDs.append((row[idColumn], qvalue))

        topNIDs = sorted(topNIDs, key=lambda x: x[1])
        input2allres[os.path.basename(pathwaysFile.name)] = [x[0] for x in topNIDs]

        for topN in args.top_n:

            localTopNIDs = [x for x in topNIDs]

            lqCount = len(localTopNIDs)
            if len(localTopNIDs) >= topN:

                topNElement = localTopNIDs[topN - 1]

                lqCount = sum([1 for x in localTopNIDs if x[1] <= topNElement[1]])

                localTopNIDs = localTopNIDs[:topN]

            #print(pathwaysFile.name, topN, len(localTopNIDs))
            foundRes[topN].append((os.path.basename(pathwaysFile.name), localTopNIDs, lqCount))


    def generate_logics(n_sets):
        """Generate intersection identifiers in binary (0010 etc)"""
        for i in range(1, 2**n_sets):
            yield bin(i)[2:].zfill(n_sets)

    def generate_petal_sets(datasets, fmt="{size}"):
        """Generate petal descriptions for venn diagram based on set sizes"""
        datasets = list(datasets)
        n_sets = len(datasets)
        dataset_union = set.union(*datasets)

        petal_sets = {}
        for logic in generate_logics(n_sets):
            included_sets = [
                datasets[i] for i in range(n_sets) if logic[i] == "1"
            ]
            excluded_sets = [
                datasets[i] for i in range(n_sets) if logic[i] == "0"
            ]
            petal_set = (
                (dataset_union & set.intersection(*included_sets)) -
                set.union(set(), *excluded_sets)
            )
            petal_sets[logic] = petal_set
        return petal_sets



    for topN in args.top_n:

        outname = args.output + "." + str(topN)
        if not outname.endswith(".png"):
            outname += ".png"

        print(outname)

        if len(foundRes[topN]) == 1:
            foundRes[topN].append(("Dummy Element",  [(None, None)], 0))

        inputSets = [set([y[0] for y in x[1]]) for x in foundRes[topN]]

        totalElements = sum([len(x) for x in inputSets])

        if totalElements == 0:
            print(outname, "no data")
            fig = plt.figure()
            plt.title("No Data to show")
            plt.savefig(outname, bbox_inches="tight")
            plt.close()

            continue

        vennLabels = venn.generate_petal_labels(inputSets)
        vennSets = generate_petal_sets(inputSets)

        for field in vennSets:

            specificField = (Counter(field)["1"] == 1)

            if not specificField:
                continue

            for specID in vennSets[field]:               
                foundIn = [ input2allres[x].index(specID) if specID in input2allres[x] else -1 for x in input2allres ]
                print(field, specID, foundIn)

        fig, ax = pwcount2function[len(foundRes[topN])](vennLabels, names=["{fn} (lq={lqc})".format(fn=x[0], lqc=x[2]) for x in foundRes[topN]])

        plt.suptitle("Overlaps for topN={} pathways (by qvalue)".format(topN))
        plt.savefig(outname, bbox_inches ="tight")


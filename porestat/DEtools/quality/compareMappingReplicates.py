
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl
import argparse
import matplotlib.patches as patches
import sys, os
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


mpl.style.use("seaborn")




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-d', '--conditions', nargs='+', type=str, action='append',  required=True, help='alignment files')
    parser.add_argument('-n', '--prefixes', nargs='+', type=str)
    parser.add_argument('-p', '--pathname', action="store_true", default=False)
    parser.add_argument('-o', '--output', type=str, required=False, help="output base")
    args = parser.parse_args()

    if not len(args.counts) == 2:
        argparse.ArgumentError("counts must have 2 files listed")

    if not len(args.conditions) == 2:
        argparse.ArgumentError("conditions must have 2 lists")

    if not len(args.prefixes) == 2:
        argparse.ArgumentError("prefixes must have 2 names")

    dfs = []
    dfHeaders = []

    countdatas = []

    for fidx, defile in enumerate(args.counts):

        print("Loading file", defile)

        indf = DataFrame.parseFromFile(defile.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })

        dfs.append(indf)

        inHeaders = indf.getHeader()
        dfHeaders.append(inHeaders)
        idColName = "Geneid" if "Geneid" in inHeaders else "id"

        # if not args.gene in inHeaders:
        #    print("Unknown gene id column", args.gene)
        #    print(inHeaders)
        #    exit(-1)

        for condition in args.conditions[fidx]:

            if not condition in inHeaders:
                print("Unknown condition", condition)
                print("Known conditions", inHeaders)
                exit(-1)

        dfData = defaultdict(lambda: dict())

        for row in indf:
            for condition in args.conditions[fidx]:
                dfData[condition][row[idColName]] = row[condition]

        countdatas.append(dfData)


    commonConditions = sorted(set.intersection(*map(set,args.conditions)))
    num_conditions = len(commonConditions)

    fig, axes = plt.subplots(num_conditions, 1, figsize=(12, 24), sharex=True, sharey=True)
    fig.suptitle("Comparison of read counts", fontsize=15)


    for i in range(0, num_conditions):

        #dfCond1 = args.conditions[0][i]
        #dfCond2 = args.conditions[1][i]
        dfCond1 = commonConditions[i]
        dfCond2 = commonConditions[i]

        cond1Values = []
        cond2Values = []

        unionGenes = set([x for x in countdatas[0][dfCond1]]+[x for x in countdatas[1][dfCond2]])

        for gene in unionGenes:

            cond1Count = countdatas[0][dfCond1].get(gene, 0) +1
            cond2Count = countdatas[1][dfCond2].get(gene, 0) +1

            cond1Values.append(cond1Count)
            cond2Values.append(cond2Count)

        print(dfCond1, dfCond2, len(cond1Values), len(cond2Values))

        axes[i].scatter(cond1Values, cond2Values, s=5)
        axes[i].set_xlabel(args.prefixes[0] + " " + dfCond1)
        axes[i].set_ylabel(args.prefixes[1] + " " + dfCond2)
        axes[i].set_yscale('log')
        axes[i].set_xscale('log')
        #axes[i].legend()

    plt.savefig(args.output + ".png", bbox_inches="tight")
    plt.close()
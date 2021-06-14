

import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl
import argparse
import math

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


#mpl.style.use("seaborn")




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-d', '--conditions', nargs='+', type=str, action='append',  required=True, help='alignment files')

    parser.add_argument('-g', '--genes', type=str, nargs='+', default=[], help="gene id column name")
    parser.add_argument('-p', '--pathname', action="store_true", default=False)

    parser.add_argument('-o', '--output', nargs='+', type=str, required=False, help="output base")

    args = parser.parse_args()

    if args.output == None:
        args.output = [counts.name for counts in args.counts]


    for fidx, defile in enumerate(args.counts):
        indf = DataFrame.parseFromFile(defile.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })

        inHeaders = indf.getHeader()

        #if not args.gene in inHeaders:
        #    print("Unknown gene id column", args.gene)
        #    print(inHeaders)
        #    exit(-1)

        for conditions in args.conditions:

            for condition in conditions:

                if not condition in inHeaders:
                    print("Unknown condition", condition)
                    print("Known conditions", inHeaders)
                    exit(-1)

        if len(args.conditions) != 2:
            print("compare 2 conditions")
            exit(-1)

        conditions1 = args.conditions[0]
        conditions2 = args.conditions[1]

        if conditions1[0].startswith("X."):
            sconditions1 = [".".join(x.split(".")[-4:-2]) for x in conditions1]
        else:
            if args.pathname:

                # just in case it's not an actual path :D
                sconditions1 = []
                for x in conditions1:
                    nx = x.split(os.path.sep)
                    sconditions1.append(nx[max([0, len(nx) - 2])])

            else:
                sconditions1 = [os.path.basename(x) for x in conditions1]


        if conditions2[0].startswith("X."):
            sconditions2 = [".".join(x.split(".")[-4:-2]) for x in conditions2]
        else:

            if args.pathname:

                # just in case it's not an actual path :D
                sconditions2 = []
                for x in conditions2:
                    nx = x.split(os.path.sep)
                    sconditions2.append(nx[max([0, len(nx) - 2])])

            else:
                sconditions2 = [os.path.basename(x) for x in conditions2]
            
        fig = plt.figure(figsize=(16, 12))

        for i in range(0, len(conditions1)):

            for j in range(0, len(conditions2)):

                cond1counts = []

                
                for ridx, row in enumerate(indf):

                    cond11count = row[conditions1[i]]
                    cond12count = row[conditions2[j]]

                    if ("Geneid" in inHeaders and row['Geneid'] in args.genes):
                        print(cond11count, cond12count)

                    if cond11count == 0:
                        cond11count = 0.1

                    if cond12count == 0:    
                        cond12count = 0.1

                    rowfc = math.log2(cond11count/cond12count)

                    if ("Geneid" in inHeaders and row['Geneid'] in args.genes):
                        print(row['Geneid'], sconditions1[i], sconditions2[j], rowfc)

                    cond1counts.append(rowfc)

                plt.hist(cond1counts, bins=len(cond1counts), label=sconditions1[i] + " vs " + sconditions2[j], density=True, cumulative=True, histtype="step")

        plt.legend(loc='upper left')

        plt.savefig(args.output[fidx] + ".interlogfc."+str(fidx)+".png", bbox_inches="tight")

        plt.close()

                        

                    
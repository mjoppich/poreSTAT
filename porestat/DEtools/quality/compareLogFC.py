

import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl
import argparse
import math

import sys, os
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


#mpl.style.use("seaborn")




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-d', '--conditions', nargs='+', type=str, action='append',  required=True, help='alignment files')
    parser.add_argument('-o', '--output', nargs='+', type=str, required=False, help="output base")
    parser.add_argument('-p', '--pathname', action="store_true", default=False)

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

        for cidx, conditions in enumerate(args.conditions):

            if conditions[0].startswith("X."):
                sconditions = [".".join(x.split(".")[-4:-2]) for x in conditions]
            else:
                if args.pathname:
                    # just in case it's not an actual path :D
                    sconditions = []
                    for x in conditions:
                        nx = x.split(os.path.sep)
                        sconditions.append(nx[max([0, len(nx)-2])])

                else:
                    sconditions = [os.path.basename(x) for x in conditions]            

            fig = plt.figure()

            for i in range(0, len(conditions)):

                for j in range(i+1, len(conditions)):

                    cond1counts = []

                    
                    for row in indf:

                        cond11count = row[conditions[i]]
                        cond12count = row[conditions[j]]

                        if cond11count == 0:
                            cond11count = 0.1

                        if cond12count == 0:    
                            cond12count = 0.1

                        rowfc = math.log2(cond11count/cond12count)
                        cond1counts.append(rowfc)

                    plt.hist(cond1counts, bins=len(cond1counts), label=sconditions[i] + " vs " + sconditions[j]  + " (n="+str(len(cond1counts))+")", density=True, cumulative=True, histtype="step")

            plt.legend(loc='upper left')

            plt.savefig(args.output[fidx] + ".logfc."+str(fidx) + "." + str(cidx) +".png", bbox_inches="tight")
            plt.close()

                        

                    
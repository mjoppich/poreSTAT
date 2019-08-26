

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
    parser.add_argument('-p', '--pathname', action="store_true", default=False)
    #parser.add_argument('-g', '--gene', type=str, required=True, help="gene id column name")
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

        for cidx, conditions in enumerate(args.conditions):

            fig, axes = plt.subplots(len(conditions), len(conditions), figsize=(8, 6), sharex=True, sharey=True)
            fig.suptitle("Comparison of read counts", fontsize=15)

            if conditions[0].startswith("X."):
                sconditions = [".".join(x.split(".")[-4:-2]) for x in conditions]
            else:

                if args.pathname:
                    sconditions = [x.split(os.path.sep)[-2] for x in conditions]
                else:
                    sconditions = [os.path.basename(x) for x in conditions]

            sconditions = ["\n".join(x.split(".")) for x in sconditions]
                
            axes[0,0].set_title(sconditions[0])

            for i in range(0, len(conditions)):

                # build a rectangle in axes coords
                left, width = .25, .5
                bottom, height = .25, .5
                right = left + width
                top = bottom + height

                # axes coordinates are 0,0 is bottom left and 1,1 is upper right
                p = patches.Rectangle(
                    (left, bottom), width, height,
                    fill=False, transform=axes[i,i].transAxes, clip_on=False
                    )

                axes[i,i].add_patch(p)

                axes[i,i].text(0.5*(left+right), 0.5*(bottom+top), sconditions[i],
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=7, color='red',
                transform=axes[i,i].transAxes)

                for j in range(i+1, len(conditions)):

                    cond1counts = []
                    cond2counts = []

                    for row in indf:

                        cond1count = row[conditions[i]]
                        cond2count = row[conditions[j]]

                        if cond1count == 0:
                            cond1count = 0.1

                        if cond2count == 0:
                            cond2count = 0.1

                        cond1counts.append(cond1count)
                        cond2counts.append(cond2count)

                    axes[i,j].scatter(cond1counts, cond2counts, label=sconditions[i] + " vs. " + sconditions[j])
                    axes[i,j].set_yscale('log')
                    axes[i,j].set_xscale('log')

                
                    if i == 0:
                        axes[i,j].set_title(sconditions[j])
                
            plt.savefig(args.output[fidx] + ".replicates."+str(fidx) + "." + str(cidx) +".png", bbox_inches ="tight")
            plt.close()
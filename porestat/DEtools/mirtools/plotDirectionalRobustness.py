import argparse
from collections import defaultdict

import matplotlib.pyplot as plt

import sys
import numpy as np
from upsetplot import from_contents, plot


def transformPVal(pval):
    return -np.log2(pval)


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-t', '--tsv', type=argparse.FileType('r'), nargs="+", required=True, help='de files')
    parser.add_argument('-p', '--prefixes', type=str, nargs='+')
    parser.add_argument('-o', '--output', type=str, required=True, help='de files')
    args = parser.parse_args()

    assert(len(args.tsv) == len(args.prefixes))

    dir2file2sigs = defaultdict(lambda: defaultdict(set))

    pvalThreshold = 0.1

    for tidx, indexFile in enumerate(args.tsv):

        for line in indexFile:
            line = line.strip().split("\t")

            elemName = line[0]
            dir = line[8]
            adjpval = line[7]

            try:
                adjpval = float(adjpval)

                if adjpval < pvalThreshold:
                    dir2file2sigs[dir][args.prefixes[tidx]].add(elemName)
            except:
                pass



    for dir in dir2file2sigs:

        allFileData = dir2file2sigs[dir]

        outname = args.output + "." + dir

        if not outname.endswith(".png"):
            outname += ".png"


        if len(allFileData) > 1:

            upIn = from_contents(allFileData)
            plot(upIn, subset_size="auto")

            plt.title("Overlap for {} elements (pval < {}".format(dir, pvalThreshold))

            plt.savefig(outname, bbox_inches="tight")
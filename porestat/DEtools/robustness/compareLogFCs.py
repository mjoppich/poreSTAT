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


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d', '--de', nargs='+', type=argparse.FileType('r'), required=True, help='DE files')
    parser.add_argument('-m', '--methods', nargs='+', type=str, required=True, help="methods to compare")

    args = parser.parse_args()
    
    for didx, detableFile in enumerate(args.de):

        indf = DataFrame.parseFromFile(detableFile.name, skipChar='#', replacements={
            "None": None,
            "": None,
            "NA": None
        })

        inHeaders = indf.getHeader()
        genesymname = None

        if "gene_symbol" in inHeaders:
            genesymname = "gene_symbol"
        elif "Geneid" in inHeaders:
            genesymname = "Geneid"
        else:
            genesymname = "id"

        print(detableFile.name)

        allMethods = args.methods

        for i in range(0, len(allMethods)):
            for j in range(i+1, len(allMethods)):

                m1 = allMethods[i]
                m2 = allMethods[j]

                m1col = m1 + "_log2FC"
                m2col = m2 + "_log2FC"

                outname = detableFile.name + ".lfc_compare.{}_{}.png".format(m1,m2)

                assert(m1col in indf.getHeader())
                assert(m2col in indf.getHeader())

                m1values = indf.getColumn(m1col)
                m2values = indf.getColumn(m2col)

                m1pvalues = []
                m2pvalues = []

                x1count = 0
                x2count = 0

                for x1,x2 in zip(m1values, m2values):

                    if not x1 is None:
                        x1count += 1
                    
                    if not x2 is None:
                        x2count += 1

                    if x1 is None or x2 is None:
                        continue

                    m1pvalues.append(float(x1))
                    m2pvalues.append(float(x2))

                minX = min(m1pvalues)
                maxX = max(m1pvalues)

                print(outname, minX, maxX)

                plt.figure()
                plt.scatter(m1pvalues, m2pvalues, color="blue", s=2, label="{} ({}) vs {} ({}) (n={} genes)".format(m1, x1count, m2, x2count, len(m1pvalues)))
                plt.plot([minX, maxX], [minX, maxX], linestyle="dashed", color="black")

                plt.xlabel(m1)
                plt.ylabel(m2)
                plt.legend(loc='lower right')
                plt.savefig(outname, bbox_inches="tight")
                plt.close()

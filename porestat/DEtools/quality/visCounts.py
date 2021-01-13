import matplotlib.pyplot as plt
from collections import defaultdict, OrderedDict
import matplotlib as mpl
import argparse

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
import pandas as pd
import seaborn as sns
import numpy as np

mpl.style.use("seaborn")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-g', '--groups', action='append', nargs='+', default=None)
    parser.add_argument('-t', '--thresholds', type=int, nargs='+', default=None)
    parser.add_argument('-o', '--output', type=str, required=False, help="output base")
    
    args = parser.parse_args()



    indf = DataFrame.parseFromFile(args.counts.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })

    dfCols = indf.getHeader()

    idColumn = None
    if "gene_symbol" in dfCols:
        idColumn = "gene_symbol"
    elif "Geneid" in dfCols:
        idColumn = "Geneid"
    elif "id" in dfCols:
        idColumn = "id"

    assert(idColumn != None)

    for grp in args.groups:
        for elem in grp:
            assert(elem in dfCols)

    alldata = OrderedDict()

    alldata["gene"] = indf.getColumn(idColumn)

    colColor = []
    colors = "rgb"

    for gidx, grp in enumerate(args.groups):

        for elem in grp:
            coldata = indf.getColumn(elem)
            alldata[elem] = coldata
            colColor.append(colors[gidx % 3])



    df = pd.DataFrame.from_dict(alldata)
    df = df.set_index("gene")
    df = (df / df.sum()) * 10000


    for t in args.thresholds:

        if t < 0:
            t = 0

        dfPlot = df[df.sum(axis=1) > t]
        dfPlot = np.log2(dfPlot+1)

        

        sns.clustermap(dfPlot, figsize=(15, 15), method="ward", col_colors = colColor, row_cluster=True)

        plt.suptitle("Heatmap of read counts (raw, library-size normalized * 10000, log2(+1), > t={} )".format(t))
        plt.savefig(args.output + ".expr_plot"+str(t)+".png")
        plt.close()


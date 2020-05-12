import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import umap
import sklearn.metrics.pairwise as pairwise
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-f', '--fc', type=argparse.FileType('r'), required=True, help='fc files')
    parser.add_argument('-o', '--output', type=str, required=True, help='output files')
    parser.add_argument('-fpkm', '--fpkm', dest='fpkm', action='store_true', default=False)
    parser.add_argument('-tpm', '--tpm', dest='tpm', action='store_true', default=False)

    parser.add_argument('-td', '--top_de', nargs="+", type=str, default=None,
                        help="instead of largest expression, top differential genes from method (must have --num)")

    parser.add_argument('-s', '--suffix', nargs='+', type=str, default=[".bam", ".sam"])
    parser.add_argument('-s2', '--samples', nargs='+', type=str, default=None)

    parser.add_argument('-n', '--num', type=int, default=-1)
    parser.add_argument('-t', '--tuple', type=int, nargs='+', default=[-3, -2])
    args = parser.parse_args()

    filename = args.fc.name

    # filename="./save/test1.tsv"
    # filename = "../../sanne_sprnaseq/data/20190703_Soehnlein/hisat2.prim.O.5.counts.fpkm.tpm.tsv"

    df = pd.read_csv(filename, skipinitialspace=True, sep='\t', comment='#')

    if args.samples != None:

        if args.tpm:
            accKeys = [x+".TPM" for x in args.samples]
        elif args.fpkm:
            accKeys = [x + ".FPKM" for x in args.samples]
        else:
            accKeys = args.samples


    else:
        suffix = args.suffix

        if args.tpm:
            suffix = [x+".TPM" for x in args.suffix]
        elif args.fpkm:
            suffix = [x + ".FPKM" for x in args.suffix]

        accKeys = [x for x in df.keys() if any([x.endswith(nsuff) for nsuff in suffix])]

    if len(accKeys) == 0:
        print("Could not find any samples")
        exit(-1)

    #for x in accKeys:
    #    print("Sample", x)

    idColName = "Geneid" if "Geneid" in df.keys() else "id"

    #print("Id Column", idColName)
    print(accKeys)

    for akey in accKeys:
        if not akey in df.keys():
            print("incorr key", akey)

    subsetDF = df[accKeys]
    subsetDF.index = df[[idColName]]

    gene2symbol = {}

    if "gene_symbol" in df.keys():

        for idx, row in df.iterrows():
            geneid = row[idColName]
            symbol = str(row["gene_symbol"])

            if len(symbol) > 0:
                gene2symbol[geneid] = symbol

        #print("Added gene2symbol", len(gene2symbol))

    if args.num != -1:

        if args.top_de != None:

            targetCols = []
            for topdeElem in args.top_de:
                targetCols += [x for x in df.keys() if x.startswith(topdeElem) and x.endswith("log2FC")]
            #print("Target cols", targetCols)

            targetColsPVal = [x.replace("log2FC", "ADJ.PVAL") for x in targetCols]
            #print("Target cols pval", targetColsPVal)

            id2val = defaultdict(lambda: 0)
            for index, row in df.iterrows():
                idVal = row[idColName]
                for cidx, col in enumerate(targetCols):

                    colVal = row[col]

                    if row[targetColsPVal[cidx]] > 0.05:
                        continue

                    if abs(colVal) > abs(id2val[idVal]):
                        id2val[idVal] = colVal

            topGenesDE = sorted([(x, id2val[x]) for x in id2val], key=lambda x: x[1], reverse=True)

            # UPREG GENES
            selGenesUp = topGenesDE[:min(len(topGenesDE), int(args.num / 2))]
            selGenesUp = [(x[0],) for x in selGenesUp if x[1] >= 0]

            selGenesDown = topGenesDE[max(len(topGenesDE) - int(args.num / 2), 0):]
            selGenesDown = [(x[0],) for x in selGenesDown if x[1] <= 0]

            #print("Upreg genes", len(selGenesUp))
            #print("Downreg genes", len(selGenesDown))

            downSubset = subsetDF.loc[subsetDF.index.isin(selGenesDown)]
            dfDown = pd.DataFrame(downSubset, columns=subsetDF.keys())

            upSubset = subsetDF.loc[subsetDF.index.isin(selGenesUp)]
            dfUp = pd.DataFrame(upSubset, columns=subsetDF.keys())

            #print(dfDown.shape)
            #print(dfUp.shape)

            subsetDF = pd.concat([dfDown, dfUp], ignore_index=False)

            #print("After up", subsetDF.shape)

            #print("after down", subsetDF.shape)
            #print(subsetDF.shape)

        else:
            subsetDF = subsetDF.nlargest(args.num, columns=accKeys)
            #print(subsetDF.shape)


    if len(gene2symbol) > 0:

        dfIndex = subsetDF.index

        symIndex = []
        for x in dfIndex:

            symIndex.append(gene2symbol.get(x[0], x[0]))

        assert (len(dfIndex) == len(symIndex))
        subsetDF.index = symIndex





    tsneDF = subsetDF
    dimNames = list(tsneDF.index)

    metric = ""

    #print("tsne shape")
    #print(tsneDF.shape)

    if tsneDF.shape[0] < 2:
        print("No data to analyse.")
        print(tsneDF.shape)
        print(tsneDF.values)
        print(tsneDF)
        exit()

    """
    HERE IS A LOG!
    """
    tsneDF = tsneDF.apply(pd.to_numeric, errors='ignore')
    tsneDF = tsneDF.replace(0, np.nan).apply(np.log10).replace(np.nan, 0)

    sns.clustermap(tsneDF, figsize=(14, 22), row_cluster=False, yticklabels=1)

    # sns.clustermap(corDF, cmap="mako", robust=True, figsize=(14,8), method='weighted', metric="correlation") #Plot the correlation as heat map
    plt.subplots_adjust(right=0.7, bottom=0.3)

    if args.top_de:
        plt.title(
            "Clustering of log-expression values from top {} down and top {} up regulated genes".format(len(selGenesDown),
                                                                                                    len(selGenesUp)))
    else:
        plt.title("Clustering of log-expression values from " + str(tsneDF.shape[0]) + " elements")

    plt.savefig(args.output + ".expr.cmap.png", bbox_inches="tight")
    plt.close()
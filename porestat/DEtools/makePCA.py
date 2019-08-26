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
    
    parser.add_argument('-cos', '--cosine', action='store_true', default=False, help="activate cosine similarity")
    parser.add_argument('-man', '--manhattan', action='store_true', default=False, help="activate manhattan similarity")
    parser.add_argument('-eucl', '--euclidean', action='store_true', default=False, help="activate manhattan similarity")

    parser.add_argument('-td', '--top_de', nargs="+", type=str, default=None, help="instead of largest expression, top differential genes from method (must have --num)")

    parser.add_argument('-s', '--suffix', type=str, default=".bam")
    parser.add_argument('-n', '--num', type=int, default=-1)
    parser.add_argument('-t', '--tuple', type=int, nargs='+', default=[-3,-2])
    args = parser.parse_args()

    filename = args.fc.name

    #filename="./save/test1.tsv"
    #filename = "../../sanne_sprnaseq/data/20190703_Soehnlein/hisat2.prim.O.5.counts.fpkm.tpm.tsv"

    df = pd.read_csv(filename, skipinitialspace=True, sep='\t', comment='#' )

    suffix = args.suffix

    if args.tpm:
        suffix += ".TPM"
    elif args.fpkm:
        suffix += ".FPKM"

    accKeys = [x for x in df.keys() if x.endswith(suffix)]

    if len(accKeys) == 0:
        print("Could not find any samples")
        exit(-1)

    idColName = "Geneid" if "Geneid" in df.keys() else "id"

    print("Id Column", idColName)

    subsetDF = df[accKeys]
    subsetDF.index = df[[idColName]]

    if args.num != -1:

        if args.top_de != None:

            targetCols = []
            for topdeElem in args.top_de:
                targetCols += [x for x in df.keys() if x.startswith(topdeElem) and x.endswith("log2FC")]
            print("Target cols", targetCols)

            targetColsPVal = [x.replace("log2FC", "ADJ.PVAL") for x in targetCols]
            print("Target cols pval", targetColsPVal)


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
            selGenesUp = topGenesDE[:min(len(topGenesDE), int(args.num/2))]
            selGenesUp = [(x[0], ) for x in selGenesUp if x[1] >= 0]

            selGenesDown = topGenesDE[max(len(topGenesDE)-int(args.num/2), 0):]
            selGenesDown = [(x[0], ) for x in selGenesDown if x[1] <= 0]

            print("Upreg genes", len(selGenesUp))
            print("Downreg genes", len(selGenesDown))
            
            odf = subsetDF

            subsetDF = odf.loc[odf.index.isin(selGenesUp)]
            subsetDF.append( odf.loc[odf.index.isin(selGenesDown)] )
            print(subsetDF.shape)

        else:
            subsetDF = subsetDF.nlargest(args.num, columns=accKeys)
            print(subsetDF.shape)

    tsneDF = subsetDF.transpose()
    dimNames = list(tsneDF.index)

    metric=""
    
    
    if args.cosine:
        cor=pairwise.cosine_similarity(tsneDF)
        cor = 1-cor
        metric="cosine"
    elif args.manhattan:
        cor=pairwise.manhattan_distances(tsneDF)
        metric="manhattan"
    elif args.euclidean:
        cor=pairwise.euclidean_distances(tsneDF)
        metric="euclidean"
    else:
        metric="correlation"
        cor=tsneDF.transpose().corr()
        cor=1-cor

    corDF = pd.DataFrame(cor)
    corDF.index = dimNames
    corDF.columns = dimNames

    linkage = hc.linkage(sp.distance.squareform(corDF, checks=False), method='weighted')
    sns.clustermap(corDF, row_linkage=linkage, col_linkage=linkage,  figsize=(14,14), )

    #sns.clustermap(corDF, cmap="mako", robust=True, figsize=(14,8), method='weighted', metric="correlation") #Plot the correlation as heat map
    plt.subplots_adjust(right=0.7, bottom=0.3)
    

    if args.top_de:
        plt.title("Clustering of expression values from top {} down and top {} up regulated genes".format(len(selGenesDown), len(selGenesUp)))
    else:
        plt.title("Clustering of expression values from "+ str(tsneDF.shape[0]) + " elements")

    plt.savefig(args.output + ".hmap.png", bbox_inches="tight")
    plt.close()


    print(dimNames)

    X_embedded = umap.UMAP(n_neighbors=3, metric=metric).fit_transform(tsneDF.values)

    plt.figure(figsize=(12,12))

    labels = dimNames

    markerTypes = ["D", "o", "^", "s", "<", ">", "1", "p", "P", "*", "H", "X", "x", "d"]

    dimTuples = []
    for x in labels:
        xs = x.split("/")

        labelTuple = []
        for pos in args.tuple:
            labelTuple.append(xs[pos])

        dimTuples.append(tuple(labelTuple))

    tupleIndices = list(set(dimTuples))
    markers = [markerTypes[tupleIndices.index(x) % len(markerTypes)] for x in dimTuples]

    for i in range(0, len(dimNames)):
        plt.scatter(X_embedded[i,0], X_embedded[i,1], label="|".join(dimTuples[i]) + " " + dimNames[i].split("/")[-1].split(".")[0], marker=markers[i])

    plt.xlabel("UMAP dim1")
    plt.ylabel("UMAP dim2")

    if args.top_de:
        plt.title("UMAP-clustering of expression values from top {} down and top {} up regulated genes for {} samples".format(len(selGenesDown), len(selGenesUp), tsneDF.shape[1]))
    else:
        plt.title("UMAP-clustering of expression values from {} elements for {} samples".format(tsneDF.shape[0], tsneDF.shape[1]))

    plt.legend()
    plt.savefig(args.output + ".umap.png", bbox_inches="tight")


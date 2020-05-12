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
    parser.add_argument('-eucl', '--euclidean', action='store_true', default=False,
                        help="activate manhattan similarity")

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

        suffix = accKeys


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

    idColName = "id"

    if "GeneName" in df.keys():
        idColName = "GeneName"
    if "Geneid" in df.keys():
        idColName = "Geneid"


    #print("Id Column", idColName)
    print(accKeys)

    subsetDF = df[accKeys]
    subsetDF.index = df[[idColName]]

    if args.top_de != None:

        targetCols = []
        for topdeElem in args.top_de:
            targetCols += [x for x in df.keys() if x.startswith(topdeElem) and x.endswith("log2FC")]
        # print("Target cols", targetCols)

        targetColsPVal = [x.replace("log2FC", "ADJ.PVAL") for x in targetCols]
        # print("Target cols pval", targetColsPVal)

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
        if args.num != -1:
            selGenesUp = topGenesDE[:min(len(topGenesDE), int(args.num / 2))]
        else:
            selGenesUp = [x for x in topGenesDE]

        selGenesUp = [(x[0],) for x in selGenesUp if x[1] >= 0]

        if args.num != -1:
            selGenesDown = topGenesDE[max(len(topGenesDE) - int(args.num / 2), 0):]
        else:
            selGenesDown = [x for x in topGenesDE]

        selGenesDown = [(x[0],) for x in selGenesDown if x[1] <= 0]

        # print("Upreg genes", len(selGenesUp))
        # print("Downreg genes", len(selGenesDown))

        dfDown = pd.DataFrame(subsetDF.loc[subsetDF.index.isin(selGenesDown)], columns=subsetDF.keys())
        dfUp = pd.DataFrame(subsetDF.loc[subsetDF.index.isin(selGenesUp)], columns=subsetDF.keys())

        # print(dfDown.shape)
        # print(dfUp.shape)

        subsetDF = pd.concat([dfDown, dfUp], ignore_index=True)

        # print("After up", subsetDF.shape)

        # print("after down", subsetDF.shape)
        # print(subsetDF.shape)

    else:

        if args.num != -1:
            subsetDF = subsetDF.nlargest(args.num, columns=accKeys)
        # print(subsetDF.shape)


    tsneDF = subsetDF.transpose()
    tsneDF = tsneDF.apply(pd.to_numeric, errors='ignore')
    
    dimNames = list(tsneDF.index)

    metric = ""

    #print("tsne shape")
    #print(tsneDF.shape)

    #print("tsne values")
    #print(tsneDF.values)


    if tsneDF.shape[1] < 2:
        print("No data to analyse.")
        print(tsneDF.shape)
        print(tsneDF.values)
        print(tsneDF)
        exit()


    if args.cosine:
        cor = pairwise.cosine_similarity(tsneDF)
        cor = 1 - cor
        metric = "cosine"
    elif args.manhattan:
        cor = pairwise.manhattan_distances(tsneDF)
        metric = "manhattan"
    elif args.euclidean:
        cor = pairwise.euclidean_distances(tsneDF)
        metric = "euclidean"
    else:
        metric = "correlation"
        cor = tsneDF.transpose().corr()

        print(cor.values)
        print(tsneDF.transpose().corr())
        print(tsneDF.transpose().dtypes)
        if np.amax(cor.values) > 1.0:
            #print("Fixing similarity")
            #print(np.amax(cor.values))

            #necessary because slightly larger values than 1.0 observed
            cor = cor / np.amax(cor.values)

        cor = 1 - cor

    corDF = pd.DataFrame(cor)
    corDF.index = dimNames
    corDF.columns = dimNames

    #print(corDF.values)


    linkage = hc.linkage(sp.distance.squareform(corDF, checks=False), method='weighted')
    sns.clustermap(corDF, row_linkage=linkage, col_linkage=linkage, figsize=(14, 14), )

    # sns.clustermap(corDF, cmap="mako", robust=True, figsize=(14,8), method='weighted', metric="correlation") #Plot the correlation as heat map
    plt.subplots_adjust(right=0.7, bottom=0.3)

    if args.top_de:
        plt.title(
            "Clustering of expression values from top {} down and top {} up regulated genes".format(len(selGenesDown),
                                                                                                    len(selGenesUp)))
    else:
        plt.title("Clustering of expression values from " + str(tsneDF.shape[0]) + " elements")

    plt.savefig(args.output + ".hmap.png", bbox_inches="tight")
    plt.close()

    #print(dimNames)

    X_embedded = umap.UMAP(n_neighbors=3, metric=metric).fit_transform(tsneDF.values)

    plt.figure(figsize=(12, 12))

    labels = dimNames

    markerTypes = ["D", "o", "^", "s", "<", ">", "1", "p", "P", "*", "H", "X", "x", "d"]

    dimTuples = []
    for x in labels:
        xs = x.split("/")

        labelTuple = []

        minTuple = min(args.tuple)
        if len(xs) < abs(minTuple):
            args.tuple = (0, )

        for pos in args.tuple:
            labelTuple.append(xs[pos])

        dimTuples.append(tuple(labelTuple))

    tupleIndices = list(set(dimTuples))
    markers = [markerTypes[tupleIndices.index(x) % len(markerTypes)] for x in dimTuples]

    for i in range(0, len(dimNames)):
        plt.scatter(X_embedded[i, 0], X_embedded[i, 1],
                    label="|".join(dimTuples[i]) + " " + dimNames[i].split("/")[-1].split(".")[0], marker=markers[i])

    plt.xlabel("UMAP dim1")
    plt.ylabel("UMAP dim2")

    if args.top_de:
        plt.title(
            "UMAP-clustering of expression values from top {} down and top {} up regulated genes for {} genes".format(
                len(selGenesDown), len(selGenesUp), tsneDF.shape[1]))
    else:
        plt.title("UMAP-clustering of expression values from {} elements for {} samples".format(tsneDF.shape[0],
                                                                                                tsneDF.shape[1]))
    plt.gca().set_aspect('equal', adjustable='box')

    plt.legend()
    plt.savefig(args.output + ".umap.png", bbox_inches="tight")

    plt.close()

    from mpl_toolkits.mplot3d import Axes3D
    pcaComponents = 2
    pca = PCA(n_components=pcaComponents)
    pca.fit(tsneDF)
    dimExplained = pca.explained_variance_ratio_

    from sklearn import decomposition

    columns = ['pca_%i' % i for i in range(pcaComponents)]
    df_pca = pd.DataFrame(pca.transform(tsneDF), columns=columns, index=tsneDF.index)
    print(df_pca.head())


    plt.figure(figsize=(12, 12))

    labels = dimNames

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
        plt.scatter(df_pca.values[i, 0], df_pca.values[i, 1],
                    label="|".join(dimTuples[i]) + " " + dimNames[i].split("/")[-1].split(".")[0], marker=markers[i])

    plt.xlabel("PCA dim1 ({:.2})".format(dimExplained[0]))
    plt.ylabel("PCA dim2 ({:.2})".format(dimExplained[1]))

    if args.top_de:
        plt.title(
            "PCA of expression values from top {} down and top {} up regulated genes for {} genes".format(
                len(selGenesDown), len(selGenesUp), tsneDF.shape[1]))
    else:
        plt.title("PCA expression values from {} elements for {} samples".format(tsneDF.shape[1],
                                                                                                tsneDF.shape[0]))
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()
    plt.savefig(args.output + ".pca.png", bbox_inches="tight")

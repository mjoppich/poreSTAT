import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import umap
import sklearn.metrics.pairwise as pairwise
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from itertools import chain

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-f', '--fc', type=argparse.FileType('r'), required=True, help='fc files')
    parser.add_argument('-o', '--output', type=str, required=True, help='output files')

    parser.add_argument('-fpkm', '--fpkm', dest='fpkm', action='store_true', default=False)
    parser.add_argument('-tpm', '--tpm', dest='tpm', action='store_true', default=False)
    parser.add_argument('-ls', '--ls', dest='ls', action='store_true', default=False)

    parser.add_argument('-cos', '--cosine', action='store_true', default=False, help="activate cosine similarity")
    parser.add_argument('-man', '--manhattan', action='store_true', default=False, help="activate manhattan similarity")
    parser.add_argument('-eucl', '--euclidean', action='store_true', default=False, help="activate manhattan similarity")

    parser.add_argument('-td', '--top_de', nargs="+", type=str, default=None, help="instead of largest expression, top differential genes from method (must have --num)")

    parser.add_argument('-s2', '--samples', nargs='+', type=str, action='append', required=True, default=None)
    parser.add_argument('-n', '--num', type=int, default=-1)
    parser.add_argument('-t', '--tuple', type=int, nargs='+', default=[-3, -2])
    args = parser.parse_args()

    filename = args.fc.name

    df = pd.read_csv(filename, skipinitialspace=True, sep='\t', comment='#')

    accKeys = []

    sample2idx = {}
    idx2samples = defaultdict(list)

    for sampleIdx, sampleGroup in enumerate(args.samples):

        if args.tpm:
            groupKeys = [x+".TPM" for x in sampleGroup]
        elif args.fpkm:
            groupKeys = [x + ".FPKM" for x in sampleGroup]
        elif args.ls:
            groupKeys = [x + ".LS" for x in sampleGroup]
        else:
            groupKeys = sampleGroup
        
        accKeys += groupKeys

        for sample in groupKeys:
            sample2idx[sample] = sampleIdx
            idx2samples[sampleIdx].append(sample)

    accKeys = [x for x in accKeys if x in df.keys()]

    if len(accKeys) == 0:
        print("Could not find any samples")
        exit(-1)

    for idx in idx2samples:
        print(idx, idx2samples[idx])


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

        else:
            print("Taking subsetDF as is.")
        # print(subsetDF.shape)

    print("subsetDF shape", subsetDF.shape)

    tsneDF = subsetDF.transpose()

    tsneDF = tsneDF.replace(to_replace=[None], value=0)
    tsneDF = tsneDF.replace(to_replace=["None"], value=0)
    tsneDF = tsneDF.apply(pd.to_numeric, errors='ignore')


    tsneDF.to_csv('/tmp/output2.tsv',sep='\t', quoting=None)
    
    dimNames = list(tsneDF.index)

    metric = ""

    if tsneDF.shape[1] < 2:
        print("No data to analyse.")
        print(tsneDF.shape)
        #print(tsneDF.values)
        #print(tsneDF)

        plt.figure()
        plt.title("No Data To Analyse (no sig genes)")
        plt.savefig(args.output + ".pca.png", bbox_inches="tight")
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

        #print(tsneDF.transpose().head())
        cor = tsneDF.transpose().corr()
        if cor.isnull().values.any():
            cor = cor.fillna(-1)
            print("Attention: NaN in correlation values!")
        

        #print(cor.values)
        #print(tsneDF.transpose().corr())
        #print(tsneDF.transpose().dtypes)
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




    markerTypes = ["D", "o", "^", "s", "<", ">", "1", "p", "P", "*", "H", "X", "x", "d"]
    
    cmap = plt.get_cmap('viridis', len(args.samples))  # matplotlib color palette name, n colors
    sampleColors = [matplotlib.colors.rgb2hex(cmap(i)[:3]) for i in range(0,len(args.samples))]


    def makePlot( reductionName, reducedEmbedding, dimNames, dimExplained):
        plt.figure(figsize=(12, 12))


        dimTuples = []
        for x in dimNames:
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
            plt.scatter(reducedEmbedding[i, 0], reducedEmbedding[i, 1], label="|".join(dimTuples[i]), marker=markers[i], color=sampleColors[sample2idx[dimNames[i]]])

        plt.xlabel("{} dim1 ({:.2})".format(reductionName, dimExplained[0]))
        plt.ylabel("{} dim2 ({:.2})".format(reductionName, dimExplained[1]))

        if args.top_de:
            plt.title(
                "{} of expression values from top {} down and top {} up regulated genes for {} genes".format(
                    reductionName, len(selGenesDown), len(selGenesUp), tsneDF.shape[1]))
        else:
            plt.title("{} expression values from {} elements for {} samples".format(reductionName, tsneDF.shape[1],
                                                                                                    tsneDF.shape[0]))
        plt.gca().set_aspect('equal', adjustable='box')
        plt.legend(bbox_to_anchor=(0, -0.2, 1, 0), loc="upper left", mode="expand", ncol=2)

        saveFile = args.output + "."+reductionName+".png"
        print(saveFile)
        plt.savefig(saveFile, bbox_inches="tight")



    umapObj = umap.UMAP(n_neighbors=3, metric=metric)
    umapEmbedding = umapObj.fit_transform(tsneDF.values)

    makePlot("UMAP", umapEmbedding, dimNames, ("-","-"))

    pcaComponents = 2
    pcaObj = PCA(n_components=pcaComponents)
    pcaEmbedding = pcaObj.fit_transform(tsneDF.values)
    dimExplained = pcaObj.explained_variance_ratio_

    makePlot("PCA", pcaEmbedding, dimNames, dimExplained)



"""
python3 /mnt/f/dev/git/poreSTAT/porestat/DEtools/de_eval/makePCA2.py --fc /mnt/t/rnaseq/Kami_aorta_colchicine/reports/mouse_colchicine_ctrl/mouse_colchicine_ctrl.combined.limmavoom.tsv --output /mnt/t/rnaseq/Kami_aorta_colchicine/reports/mouse_colchicine_ctrl/mouse_colchicine_ctrl.combined.limmavoom.tsv.LS.all_expr.mpca --num -1 --ls --samples readsexon_colchicine_1 readsexon_colchicine_2 readsexon_colchicine_3 readsexon_colchicine_6 readsexon_colchicine_7 readsexon_colchicine_10 readsexon_colchicine_13 readsexon_colchicine_14 readsexon_colchicine_17 readsexon_colchicine_18 readsexon_colchicine_21 readsexon_colchicine_22 readsexon_colchicine_25 readsexon_colchicine_26 --samples umiexon_colchicine_1 umiexon_colchicine_2 umiexon_colchicine_3 umiexon_colchicine_6 umiexon_colchicine_7 umiexon_colchicine_10 umiexon_colchicine_13 umiexon_colchicine_14 umiexon_colchicine_17 umiexon_colchicine_18 umiexon_colchicine_21 umiexon_colchicine_22 umiexon_colchicine_25 umiexon_colchicine_26 --samples readsexon_control_4 readsexon_control_5 readsexon_control_8 readsexon_control_9 readsexon_control_11 readsexon_control_12 readsexon_control_15 readsexon_control_16 readsexon_control_19 readsexon_control_20 readsexon_control_23 readsexon_control_24 readsexon_control_27 readsexon_control_28 --samples umiexon_control_4 umiexon_control_5 umiexon_control_8 umiexon_control_9 umiexon_control_11 umiexon_control_12 umiexon_control_15 umiexon_control_16 umiexon_control_19 umiexon_control_20 umiexon_control_23 umiexon_control_24 umiexon_control_27 umiexon_control_28

"""
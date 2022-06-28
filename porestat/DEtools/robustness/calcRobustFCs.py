import math

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")


import matplotlib as mpl
mpl.use('Agg')


import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd

import argparse
from upsetplot import from_contents,plot
import numpy as np
from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from adjustText import adjust_text

mpl.style.use("seaborn")



def getColorMap(colormap="Viridis"):
    cmap = plt.cm.get_cmap(colormap)
    return cmap

def getColor(colormap="Viridis", value=0.5):
    cmap = getColorMap(colormap=colormap)
    return cmap(value)

def getColorLin( min, max, val, colormap="viridis"):
    value = val / (max-min)
    return getColor(colormap=colormap, value=value)

def getColorVector(cnt, colormap="Viridis"):

    colors = []
    for i in range(0, cnt):
        color = getColorLin(0, cnt, i)
        colors.append(color)
        
    return colors



def plot_volcano(FcPvalGene, title, outfile, minpval, minfc, showGeneCount=30):



    color1 = "#883656"  #"#BA507A"
    color1_nosig = "#BA507A"
    color1_nosig_less = "#d087a4"

    color2 = "#4d6841"
    color2_nosig = "#70975E"
    color2_nosig_less = "#99b78b"



    colors = {"down": (color1, color1_nosig, color1_nosig_less), "up": (color2, color2_nosig,color2_nosig_less)}
    
    with plt.style.context("default"):
        plt.figure(figsize=(16,10))

        FcPvalGene = sorted(FcPvalGene, key=lambda x: x[1])

        xydots = [(x[0], -np.log10(x[1])) for x in FcPvalGene]
        maxally = max([x[1] for x in xydots if not np.isinf(x[1])])
        xydots = [(x, y if y <= maxally else maxally) for x,y in xydots]
        dotgene = [x[2] for x in FcPvalGene]

        pvalThresh = -np.log10( minpval )

        showGeneCount_pos = showGeneCount
        showGeneCount_neg = showGeneCount

        showGenes = []
        for x in FcPvalGene:

            gene = x[2]
            geneFC = x[0]

            if geneFC < 0 and showGeneCount_neg > 0:
                showGenes.append(gene)
                showGeneCount_neg -= 1

            if geneFC > 0 and showGeneCount_pos > 0:
                showGenes.append(gene)
                showGeneCount_pos -= 1


        texts = []
        sel_down_xy = []
        nosig_down_xy = []
        nosigless_down_xy = []

        sel_up_xy = []
        nosig_up_xy = []
        nosigless_up_xy = []

        upregCount = 0
        downregCount = 0
        upregSigCount = 0
        downregSigCount = 0
        unregCount = 0

        for gi, (x,y) in enumerate(xydots):

            if x < 0:
                if y >= pvalThresh and abs(x) >= minfc:
                    downregSigCount += 1
                else:
                    downregCount += 1
            elif x > 0:
                if y >= pvalThresh and abs(x) >= minfc:
                    upregSigCount += 1
                else:
                    upregCount += 1
            elif x == 0:
                unregCount += 1

            if dotgene[gi] in showGenes:

                if x < 0:
                    sel_down_xy.append((x,y))
                else:
                    sel_up_xy.append((x,y))


                texts.append(plt.text(x * (1 + 0.01), y * (1 + 0.01) , dotgene[gi], fontsize=12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.7)))

            else:

                if x < 0:

                    if y >= pvalThresh and abs(x) >= minfc:
                        nosig_down_xy.append((x,y))
                    else:
                        nosigless_down_xy.append((x,y))
                else:
                    if y >= pvalThresh and abs(x) >= minfc:
                        nosig_up_xy.append((x,y))
                    else:
                        nosigless_up_xy.append((x,y))


        #print(len(sel_xy), "of", len(genes))
        plt.plot([x[0] for x in nosigless_up_xy], [x[1] for x in nosigless_up_xy], '.', color=colors["up"][2])
        plt.plot([x[0] for x in nosigless_down_xy], [x[1] for x in nosigless_down_xy], '.', color=colors["down"][2])

        plt.plot([x[0] for x in nosig_up_xy], [x[1] for x in nosig_up_xy], 'o', color=colors["up"][1])
        plt.plot([x[0] for x in nosig_down_xy], [x[1] for x in nosig_down_xy], 'o', color=colors["down"][1])

        plt.plot([x[0] for x in sel_up_xy], [x[1] for x in sel_up_xy], 'o', color=colors["up"][0])
        plt.plot([x[0] for x in sel_down_xy], [x[1] for x in sel_down_xy], 'o', color=colors["down"][0])


        plt.hlines(y=pvalThresh, xmin=plt.xlim()[0], xmax=-minfc, linestyle="dotted")
        plt.hlines(y=pvalThresh, xmin=minfc, xmax=plt.xlim()[1], linestyle="dotted")

        yMaxLim = plt.ylim()[1]

        plt.vlines(x=-minfc, ymin=pvalThresh, ymax=yMaxLim, linestyle="dotted")
        plt.vlines(x=minfc, ymin=pvalThresh, ymax=yMaxLim, linestyle="dotted")

        adjust_text(texts, force_points=0.2, force_text=0.2, expand_points=(2, 2), expand_text=(1, 1), arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
        #        texts.append(plt.text(x * (1 + 0.01), y * (1 + 0.01) , dotgene[gi], fontsize=12))

        plt.title(title, fontsize = 40)
        plt.xlabel("logFC", fontsize = 32)
        plt.ylabel("Neg. Log. Adj. P-Value", fontsize = 32)
        plt.xticks(fontsize=14)

        infoText = "Total Genes: {totalCount}; Down-Reg. (sig.): {downregCount} ({downregSigCount}); Up-Reg. (sig.): {upregCount} ({upregSigCount}); Un-Reg.: {unregCount}".format(
            totalCount=upregCount+downregCount+upregSigCount+downregSigCount+unregCount,
            upregCount=upregCount, upregSigCount=upregSigCount,
            downregCount=downregCount, downregSigCount=downregSigCount,
            unregCount=unregCount
        )
        plt.figtext(0.5, 0.01, infoText, wrap=True, horizontalalignment='center', fontsize=14)

        print(outfile)
        plt.savefig(outfile, bbox_inches="tight")



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d', '--de', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-m', '--methods', nargs='+', type=str, required=True, help='alignment files')
    parser.add_argument('-o', '--output', nargs='+', type=argparse.FileType('w'), required=True, help='alignment files')
    parser.add_argument('-p', '--print', action="store_true", default=False)


    #parser.add_argument('-c', '--cutoff', type=float, help='alignment files', default=0.05)
    parser.add_argument('-minfc', '--min-foldchange', type=float, default=1.0, required=False)
    parser.add_argument('-minpval', '--min-pvalue', type=float, default=0.05, required=False)


    args = parser.parse_args()


    for fidx, defile in enumerate(args.de):
        indf = DataFrame.parseFromFile(defile.name)

        availMethods = set()

        headername2idx = {}

        indfHeader = indf.getHeader()
        genesymname = None

        if "gene_symbol" in indfHeader:
            genesymname = "gene_symbol"
        elif "Geneid" in indfHeader:
            genesymname = "Geneid"
        else:
            genesymname = "id"

        col2idx = {x: indf.getColumnIndex(x) for x in indfHeader}
        for x in col2idx:
            print(x, col2idx[x])

        methodCol2Idx = defaultdict(list)

        for hidx, headername in enumerate(indfHeader):

            hsplit = headername.split("_")

            #if len(hsplit) != 2:
            #    continue

            if hsplit[-1] in ["RAW.PVAL", "ADJ.PVAL", "log2FC"]:
                availMethods.add(hsplit[0])
                methodCol2Idx[(hsplit[0], hsplit[-1])].append(hidx)

        for method in args.methods:
            if not method in availMethods:
                print("Invalid Method:", method)
                print("Available Methods:", availMethods)
                exit(-1)

        rowFCidx = indf.addColumn("ROB_log2FC")
        robRawPidx = indf.addColumn("ROB_RAW.PVAL")
        robAdjPidx = indf.addColumn("ROB_ADJ.PVAL")

        signCount = 0
        method2genes = defaultdict(set)

        def makeRobFC(line):

            global signCount
            global method2genes

            methodsAdjPvals = []
            methodsRawPvals = []
            methodsLog2FCs = []

            robustAdjPval = 1
            robustRawPval = 1
            robustLog2FC = 0

            

            for method in availMethods:

                adjPVals = methodCol2Idx[(method, "ADJ.PVAL")]
                rawPVals = methodCol2Idx[(method, "RAW.PVAL")]
                logFCs = methodCol2Idx[(method, "log2FC")]

                for i in range(0, len(adjPVals)):

                    adjPval = line[adjPVals[i]]

                    if i < len(rawPVals):
                        rawPval = line[rawPVals[i]]
                    else:
                        rawPval = None

                    logFC = line[logFCs[i]]

                    if adjPval == None or adjPval == "None" or adjPval == "NA":
                        adjPval = 1.0

                    if rawPval == None or rawPval == "None" or rawPval == "NA":
                        rawPval = 1.0

                    if logFC == None or logFC == "None" or logFC == "NA":
                        logFC = 0.0
                    
                    
                    if float(adjPval) < args.min_pvalue and abs(float(logFC)) > args.min_foldchange:
                        method2genes[method].add( line[col2idx[genesymname]] )

                    if not method in args.methods:
                        continue

                    methodsAdjPvals.append(float(adjPval))
                    methodsRawPvals.append(float(rawPval))
                    methodsLog2FCs.append(float(logFC))

            allUp = all([x >= 0 for x in methodsLog2FCs])
            allDown = all([x <= 0 for x in methodsLog2FCs])

            allSameDirection = allUp or allDown

            robustAdjPval = max(methodsAdjPvals)
            robustRawPval = max(methodsRawPvals)

            if allSameDirection:

                if allUp:
                    robustLog2FC = min(methodsLog2FCs)
                else:
                    robustLog2FC = max(methodsLog2FCs)


            else:
                robustLog2FC = 0
                robustRawPval = 1
                robustAdjPval = 1

                    #print("Inconsistent lfc", line[col2idx["gene_symbol"]], methodsLog2FCs)

            #if line[0] == 'ENSMUSG00000041773':
            #    print(methodsAdjPvals, methodsRawPvals, methodsLog2FCs)
            #    print(robustAdjPval, robustRawPval, robustLog2FC, allSameDirection)

            if robustAdjPval < args.min_pvalue  and abs(robustLog2FC) > args.min_foldchange:
                signCount += 1

                if args.print:
                    print(line[col2idx["id"]], line[col2idx["gene_symbol"]], methodsAdjPvals, methodsRawPvals, methodsLog2FCs)

            assert (robustLog2FC != None)
            assert (robustRawPval != None)
            assert (robustAdjPval != None)

            line[rowFCidx] = robustLog2FC
            line[robRawPidx] = robustRawPval
            line[robAdjPidx] = robustAdjPval

            return line


        indf.applyToRow(makeRobFC)

        fcPvals = []
        fcPvalsSig = []

        fcPvalGene = []

        for row in indf:

            rowFC = row["ROB_log2FC"]
            rowPV = row["ROB_ADJ.PVAL"]
            rowGene = row[genesymname]

            if rowFC == None or rowPV == None:
                continue

            if rowPV == 0:
                logRowPV = -1000
            else:
                try:
                    logRowPV = -math.log10(rowPV)
                except:
                    print(row.to_tuple())
                    print("logfc", row["ROB_log2FC"], "adj.pval", row["ROB_ADJ.PVAL"])
                    exit(-1)

            fcPvalGene.append((rowFC, rowPV, rowGene))

            if rowPV < args.min_pvalue and abs(rowFC) >= args.min_foldchange:
                fcPvalsSig.append((rowFC, logRowPV))
            else:
                fcPvals.append((rowFC, logRowPV))

        #plt.scatter([x[0] for x in fcPvals], [x[1] for x in fcPvals], color="blue", label="Gene (n={})".format(len(fcPvals)))
        #plt.scatter([x[0] for x in fcPvalsSig], [x[1] for x in fcPvalsSig], color="red", label="Gene pVal < 0.05 &\n logFC >= 1 (n={})".format(len(fcPvalsSig)))
        #plt.legend(loc='lower left')
        #plt.savefig(args.output[fidx].name + ".rob.volcano.png", bbox_inches="tight")

        volcanoOutfile = args.output[fidx].name + ".rob.volcano.png"
        print(volcanoOutfile)
        plot_volcano(fcPvalGene, ";".join([str(x) for x in args.methods]), volcanoOutfile, args.min_pvalue, args.min_foldchange)

        print("Writing out DF", args.output[fidx].name)
        indf.export(args.output[fidx].name)

        # explicitly sets an empty set
        for method in availMethods:
            if not method in method2genes:
                method2genes[method] = set()

        totalGenes = 0
        for x in method2genes:
            print(x, len(method2genes[x]))
            totalGenes += len(method2genes[x])

        outfilename = args.output[fidx].name + ".rob.upset.png"
        print(outfilename)
        if totalGenes == 0:
            fig = plt.figure()
            plt.title("No Data to show")
            plt.savefig(outfilename, bbox_inches="tight")
            plt.close()
        
        else:

            if len(method2genes) == 1:
                method2genes["N/A"] = {}

            plt.title("Overlap of called genes per method")
            upIn = from_contents(method2genes)
            #print( upIn.index )
            #print( upIn.index.levels )

            plot(upIn, subset_size="auto") 
            plt.savefig(outfilename, bbox_inches="tight")


    

    
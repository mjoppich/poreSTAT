import matplotlib as mpl
mpl.use('Agg')


import matplotlib.pyplot as plt
from collections import defaultdict


import argparse
from upsetplot import from_contents,plot
import sys
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


mpl.style.use("seaborn")

def makeplot(plotdata, filename, stats):

    maxlen = 0

    commonGenes = defaultdict(set)

    for sample in plotdata:
        genecount = plotdata[sample]

        for x in genecount:
            commonGenes[sample].add(x[0])


    allGenes = set()
    for s in commonGenes:
        allGenes = allGenes.union(commonGenes.get(s, set()))

    for s in commonGenes:
        allGenes = allGenes.intersection(commonGenes.get(s, set()))

    print("Common Genes in all Samples", len(allGenes))

    if len(allGenes) < 100:
        for x in allGenes:
            print(x)

    colors = getColorVector(len(plotdata))

    maxlen = 0
    for sample in plotdata:
        maxlen = max([maxlen, len(plotdata[sample])])

    print()
    print("All Counts")
    print()

    plt.figure()

    for idx, sample in enumerate(plotdata):
        counts = [x[1] for x in plotdata[sample]]
        acounts = len([x for x in counts if x >= 1+pseudoCount])
        print(sample, len(counts), acounts)
        stats[sample]["genes_with_count"] = len(counts)
        plt.hist(counts, int(maxlen), color=colors[idx], normed=False,cumulative =True, label=sample + " ("+str(len(counts))+", "+str(acounts)+")", histtype="step" )


    plt.legend()
    plt.xscale("log")
    plt.title("Histogram of read counts with any count")
    plt.savefig(filename + ".countplot.png")
    #plt.show()

    plt.close()

    minCount = 1

    print()
    print("All Counts pseudocount+"+str(minCount))
    print()

    plt.figure()

    for idx, sample in enumerate(plotdata):
        counts = [x[1] for x in plotdata[sample] if x[1] >= minCount+pseudoCount]
        print(sample, len(counts), len([x for x in counts if x >= minCount+pseudoCount]))

        stats[sample]["genes_with_mincount_" + str(minCount)] = len(counts)

        plt.hist(counts, int(maxlen), color=colors[idx], normed=False,cumulative =True, label=sample + " ("+str(len(counts))+")", histtype="step" )

    plt.legend()
    plt.xscale("log")
    plt.title("Histogram of read counts with count >=" + str(minCount))
    plt.savefig(filename + ".countplot"+str(minCount)+".png")
    plt.close()
    minCount = 20

    print()
    print("All Counts pseudocount+"+str(minCount))
    print()

    plt.figure()

    for idx, sample in enumerate(plotdata):
        counts = [x[1] for x in plotdata[sample] if x[1] >= minCount+pseudoCount]
        print(sample, len(counts), len([x for x in counts if x >= minCount+pseudoCount]))
        stats[sample]["genes_with_mincount_" + str(minCount)] = len(counts)
        plt.hist(counts, int(maxlen), color=colors[idx], normed=False,cumulative =True, label=sample + " ("+str(len(counts))+")", histtype="step" )

    plt.legend()
    plt.xscale("log")
    plt.title("Histogram of read counts with count >=" + str(minCount))
    plt.savefig(filename + ".countplot"+str(minCount)+".png")
    plt.close()

    return stats

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




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d', '--de', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-m', '--methods', nargs='+', type=str, required=True, help='alignment files')
    parser.add_argument('-o', '--output', nargs='+', type=argparse.FileType('w'), required=True, help='alignment files')
    parser.add_argument('-p', '--print', action="store_true", default=False)
    parser.add_argument('-c', '--cutoff', type=float, help='alignment files', default=0.05)


    args = parser.parse_args()


    for fidx, defile in enumerate(args.de):
        indf = DataFrame.parseFromFile(defile.name)

        availMethods = set()

        headername2idx = {}

        indfHeader = indf.getHeader()

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

        print(availMethods)
        print(args.methods)

        for method in args.methods:
            if not method in availMethods:
                print("Invalid Method:", method)
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
                    rawPval = line[rawPVals[i]]
                    logFC = line[logFCs[i]]

                    if adjPval == None or adjPval == "None" or adjPval == "NA":
                        adjPval = 1.0

                    if rawPval == None or rawPval == "None" or rawPval == "NA":
                        rawPval = 1.0

                    if logFC == None or logFC == "None" or logFC == "NA":
                        logFC = 0.0
                    
                    
                    if float(adjPval) < args.cutoff:
                        method2genes[method].add(line[0])

                    if not method in args.methods:
                        continue

                    methodsAdjPvals.append(float(adjPval))
                    methodsRawPvals.append(float(rawPval))
                    methodsLog2FCs.append(float(logFC))
        
            allSameDirection = all([x >= 0 for x in methodsLog2FCs]) or all([x <= 0 for x in methodsLog2FCs]) 

            if allSameDirection:

                robustAdjPval = max(methodsAdjPvals)
                robustRawPval = max(methodsRawPvals)

                if methodsLog2FCs[0] >= 0:
                    robustLog2FC = min(methodsLog2FCs)
                else:
                    robustLog2FC = max(methodsLog2FCs)

            #if line[0] == 'ENSMUSG00000041773':
            #    print(methodsAdjPvals, methodsRawPvals, methodsLog2FCs)
            #    print(robustAdjPval, robustRawPval, robustLog2FC, allSameDirection)

            if robustAdjPval < args.cutoff:
                signCount += 1

                if args.print:
                    print(line[col2idx["id"]], line[col2idx["gene_symbol"]], methodsAdjPvals, methodsRawPvals, methodsLog2FCs)

            line[rowFCidx] = robustLog2FC
            line[robRawPidx] = robustRawPval
            line[robAdjPidx] = robustAdjPval

            return line


        indf.applyToRow(makeRobFC)

        print(signCount)

        indf.export(args.output[fidx].name)

        for method in availMethods:
            if not method in method2genes:
                method2genes[method] = set()

        for x in method2genes:
            print(x, len(method2genes[x]))

        upIn = from_contents(method2genes)
        plot(upIn, subset_size="auto") 

        plt.savefig(args.output[fidx].name + ".upset.png", bbox_inches="tight")


    

    
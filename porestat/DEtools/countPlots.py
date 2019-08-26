import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl
import argparse

import sys
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


mpl.style.use("seaborn")

def makeplot(plotdata, filename, stats, outname):

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
    plt.savefig(outname + ".countplot.png")
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
    plt.savefig(outname + ".countplot"+str(minCount)+".png")
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
    plt.savefig(outname + ".countplot"+str(minCount)+".png")
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
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-f', '--fname',action='store_true', default=False, help='alignment files')
    parser.add_argument('-g', '--groups', action='append', nargs='+', default=None)

    parser.add_argument('-o', '--output', nargs='+', type=str, required=False, help="output base")
    
    args = parser.parse_args()

    if args.output == None:
        args.output = [counts.name for counts in args.counts]

    pseudoCount = 1

    allSamples = []
    for x in args.groups:
        for y in x:
            if not y in allSamples:
                allSamples.append(y)

    print(allSamples)

    for fidx, defile in enumerate(args.counts):
        indf = DataFrame.parseFromFile(defile.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })

        dfCols = indf.getHeader()

        allDerivedSamples = []
        replaceSamples = {}
        for sample in allSamples:

            if sample in dfCols:
                replaceSamples[sample] = sample
                allDerivedSamples.append(sample)

            else:
                sampleAdded = False
                for colName in dfCols:
                    if colName.endswith(sample[1:]):
                        sampleAdded = True
                        replaceSamples[sample] = colName
                        allDerivedSamples.append(colName)
                        break
                    elif colName.replace("//", "/").endswith(sample.replace("//", "/")[1:]):
                        sampleAdded = True
                        replaceSamples[sample] = colName
                        allDerivedSamples.append(colName)
                        break

                if not sampleAdded:
                    print("Could not assign sample", sample)
                    print(dfCols)

        allSamples = allDerivedSamples

        for sample in replaceSamples:
            print(sample, "-->", replaceSamples[sample])

        dfGroups = args.groups
        newgroups = []
        for group in dfGroups:
            ngroup = []
            for gelem in group:
                ngroup.append(replaceSamples[gelem])
            
            newgroups.append(ngroup)

        dfGroups = newgroups

        sample2genecount = defaultdict(list)

        idxes = None

        for row in indf:
            for sname in allSamples:
                rcount = float(row[sname])+pseudoCount
                sample2genecount[sname].append( (row[sname], rcount) )

        print("Total Sample Counts")

        sample2stats = {}

        for sample in sample2genecount:

            totalCount = sum([x[1]-pseudoCount for x in sample2genecount[sample]])
            print(sample, totalCount)

            sample2stats[sample] = {"sample": sample, "totalCount": totalCount}

        sample2stats = makeplot(sample2genecount, defile.name, sample2stats, args.output[fidx])

        columns = list()
        for sample in sample2stats:
            for x in sample2stats[sample]:
                if not x in columns:
                    columns.append(x)

        outdf = DataFrame()
        outdf.addColumns(columns)

        for sample in sample2stats:

            dr = DataRow.fromDict(sample2stats[sample])
            outdf.addRow(dr)

        print(outdf)

        if dfGroups != None:

            allGenes = set()
            for sample in sample2genecount:
                for x in sample2genecount[sample]:
                    allGenes.add(x[0])

            for group in dfGroups:
                for gm in group:
                    if not gm in sample2genecount:
                        print([x for x in sample2genecount])
                        print("Invalid group member name:", gm)
                        exit(-1)

            for minCount in [1, 10, 20, 100]:
                for group in dfGroups:

                    groupGenes = set([x for x in allGenes])
                    for gm in group:

                        groupGenes = groupGenes.intersection(set([x[0] for x in sample2genecount[gm] if x[1] >= minCount + pseudoCount]))

                    print(";".join(group), minCount, len(allGenes), len(groupGenes))

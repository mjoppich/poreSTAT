

import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import matplotlib as mpl
import argparse
import math

import sys, os
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


#mpl.style.use("seaborn")




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-d', '--conditions', nargs='+', type=str, action='append',  required=True, help='alignment files')
    parser.add_argument('-l', '--last', action="store_true", default=False)
    parser.add_argument('-i', '--ignoreMissing', action="store_true", default=False)
    parser.add_argument('-o', '--output', nargs='+', type=str, required=False, help="output base")
    
    args = parser.parse_args()

    if args.output == None:
        args.output = [counts.name for counts in args.counts]

    allSamples = []
    for x in args.conditions:
        for y in x:
            if not y in allSamples:
                allSamples.append(y)

    for fidx, defile in enumerate(args.counts):
        indf = DataFrame.parseFromFile(defile.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })

        inHeaders = indf.getHeader()

        newconds = args.conditions


        dfCols = indf.getHeader()

        allDerivedSamples = []
        replaceSamples = {}
        for sample in allSamples:

            if sample in inHeaders:
                replaceSamples[sample] = sample
                allDerivedSamples.append(sample)

            else:
                sampleAdded = False
                for colName in inHeaders:
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
                    print(inHeaders)

                    if not args.ignoreMissing:
                        exit(-1)

        allSamples = allDerivedSamples

        for sample in replaceSamples:
            print(sample, "-->", replaceSamples[sample])

        dfGroups = args.conditions
        newgroups = []
        for group in dfGroups:
            ngroup = []
            for gelem in group:
                if gelem in replaceSamples: # may only not be included, if not added -> error before
                    ngroup.append(replaceSamples[gelem])
            
            newgroups.append(ngroup)

        newconds = newgroups








        lastconds = []

        for conditions in newconds:

            lcond = []

            for condition in conditions:

                if args.last:

                    for x in inHeaders:

                        if condition.split("/")[-1] in x:
                            lcond.append(x)
                            break

                elif not condition in inHeaders:
                    print("Unknown condition", condition)
                    print("Known conditions", inHeaders)
                    exit(-1)

            if args.last:
                lastconds.append(lcond)

        if not args.last:
            lastconds = newconds

        condition2gc = defaultdict(list)
        total2condition = Counter()

        allconditions = set()
        for x in lastconds:
            for y in x:
                allconditions.add(y)


        for ridx, row in enumerate(indf):
            for condition in allconditions:

                counts = float(row[condition])

                total2condition[condition] += counts

                if "gene_symbol" in inHeaders:
                    condition2gc[ condition ].append( (row["gene_symbol"], counts) )
                elif "Geneid" in inHeaders:
                    condition2gc[ condition ].append( (row["Geneid"], counts) )
                elif "id" in inHeaders:
                    condition2gc[ condition ].append( (row["id"], counts) )
                
                else:
                    condition2gc[ condition ].append( (ridx, counts) )



        for cidx, conditions in enumerate(lastconds):

            fig, axes = plt.subplots(len(conditions), 1, figsize=(14, 30))

            #print(conditions)
            for condition in conditions:
                print(condition, total2condition[condition])

            for ccidx, condition in enumerate(conditions):
                sgenes = sorted(condition2gc[condition], key=lambda x: x[1],reverse=True) 
                
                printGenes = None

                if len(sgenes) < 50:
                    printGenes = sgenes

                else:
                    printGenes = sgenes[0:50]

                #print(len(printGenes))
                
                labels = [x[0] for x in printGenes]
                counts = [x[1]/total2condition[condition] for x in printGenes]

                axes[ccidx].bar(range(0, len(labels)), counts, align="center")
                axes[ccidx].set_xlabel("Gene Name")
                axes[ccidx].set_ylabel("Count Frequency")
                axes[ccidx].set_xticklabels(labels)
                axes[ccidx].set_title(condition)
                axes[ccidx].set_xticks([i for i in range(0, len(printGenes))])

                for tick in axes[ccidx].get_xticklabels():
                    tick.set_rotation(90)
        
            plt.savefig(args.output[fidx] + ".cpergenes."+str(fidx)+"." + str(cidx)+ ".png", bbox_inches ="tight")
            plt.close()

        

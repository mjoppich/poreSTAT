

import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl
import argparse
import math

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


#mpl.style.use("seaborn")




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-d', '--de', nargs='+', type=argparse.FileType('r'),  required=True, help='alignment files')
    parser.add_argument('-g', '--conditions', nargs='+', type=str, action='append',  required=True, help='alignment files')

    parser.add_argument('-p', '--pval', type=float, default=0.05)
    parser.add_argument('-t', '--tools', nargs='+')

    parser.add_argument('-o', '--output', type=argparse.FileType("w"), required=True)

    #parser.add_argument('-g', '--gene', type=str, required=True, help="gene id column name")


    args = parser.parse_args()


    for fidx, defile in enumerate(args.counts):
        indf = DataFrame.parseFromFile(defile.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })

        inHeaders = indf.getHeader()

        #if not args.gene in inHeaders:
        #    print("Unknown gene id column", args.gene)
        #    print(inHeaders)
        #    exit(-1)

        allconditions = []

        for conditions in args.conditions:

            for condition in conditions:

                if not condition in inHeaders:
                    print("Unknown condition", condition)
                    print("Known conditions", inHeaders)
                    exit(-1)
                else:
                    if not condition in allconditions:
                        allconditions.append(condition)

        dfGeneName = "id"

        if not "id" in indf.getHeader():
            dfGeneName = "gene"

        gene2samplecounts = defaultdict(lambda: dict())
        for row in indf:

            for condition in allconditions:
                gene2samplecounts[row[dfGeneName]][condition] = float(row[condition])


        dedf = indf = DataFrame.parseFromFile(args.de[fidx].name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })

        method1Only = set()
        method2Only = set()

        deGeneName = "id"

        if not "id" in dedf.getHeader():
            deGeneName = "gene"
        
        print("genename", args.tools[0], args.tools[1], args.conditions[0], args.conditions[1], sep="\t", file=args.output)
        
        for row in dedf:

            method1Pval = row[args.tools[0] + "_ADJ.PVAL"]
            method2Pval = row[args.tools[1] + "_ADJ.PVAL"]

            if method1Pval == None or method1Pval == "None" or method1Pval == "NA":
                method1Pval = 1.0
            if method2Pval == None or method2Pval == "None" or method2Pval == "NA":
                method2Pval = 1.0

            method1Pval = float(method1Pval)
            method2Pval = float(method2Pval)

            geneName = row[deGeneName]

            if (method1Pval < args.pval and not method2Pval < args.pval):
                if not geneName in gene2samplecounts:
                    #print("Unknown gene", geneName)

                    print(geneName, True, False, "", "", sep="\t", file=args.output)
                    continue

                method1Only.add(geneName)

            elif (method2Pval < args.pval and not method1Pval < args.pval):
                if not geneName in gene2samplecounts:
                    # print("Unknown gene", geneName)

                    print(geneName, False, True, "", "", sep="\t", file=args.output)
                    continue

                method2Only.add(geneName)

        for gene in method1Only.union(method2Only):

            method1counts = []
            method2counts = []

            for cidx, conditions in enumerate(args.conditions):

                if conditions[0].startswith("X."):
                    sconditions = [".".join(x.split(".")[-4:-2]) for x in conditions]
                else:
                    sconditions = [os.path.basename(x) for x in conditions]

                
                for condition in conditions:

                    if cidx == 0:
                        method1counts.append(gene2samplecounts[gene][condition])
                    elif cidx == 1:
                        method2counts.append(gene2samplecounts[gene][condition])

            
            print(gene, gene in method1Only, gene in method2Only, method1counts, method2counts, sep="\t", file=args.output)
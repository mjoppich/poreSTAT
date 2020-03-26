import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl
import argparse
import matplotlib.patches as patches
import sys, os

from statsmodels.stats.multitest import multipletests

sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d', '--de', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-s', '--sets', type=argparse.FileType('r'), required=True, help='set file')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True, help='set file')

    parser.add_argument('--pval', type=str, required=False, default='ROB_ADJ.PVAL')
    parser.add_argument('--logfc', type=str, required=False, default='ROB_log2FC')
    parser.add_argument('--gene_symbol', type=str, required=False, default='gene_symbol')

    args = parser.parse_args()

    sym2approvSym = {}

    with open("/mnt/d/dev/data/genomes/hgnc2sym2ens2uniprot") as fin:

        for line in fin:
            line = line.strip().split("\t")

            sym = line[0]
            approvSym = line[3]

            sym2approvSym[sym] = approvSym

    availSets = {}

    setDF = DataFrame.parseFromFile(args.sets, skipChar='#', replacements={
        "None": None,
        "": None,
        "NA": None
    })

    allSetGenes = set()
    for row in setDF:
        availSets[row['set_id']] = (row["set_descr"], set([sym2approvSym[x.strip()] for x in row["genes"].split(";") if x.strip() in sym2approvSym]))
        allSetGenes = allSetGenes.union(availSets[row['set_id']][1])

    print("Got", len(availSets), "sets with a total of", len(allSetGenes), "genes")

    indf = DataFrame.parseFromFile(args.de.name, skipChar='#', replacements={
        "None": None,
        "": None,
        "NA": None
    })

    inHeaders = indf.getHeader()

    outdf = DataFrame()
    outdf.addColumns(['elem_id', 'population_size', 'success_population', 'sample_size', 'success_samples', 'sample_success_fraction', 'pval', 'adj_pval', 'direction', 'genes'])

    for direction in ["UP", "DOWN", "ANY"]:

        significantGenes = set()
        measuredGenes = set()

        for row in indf:

            geneID = row[args.gene_symbol]

            if geneID in sym2approvSym:
                geneID = sym2approvSym[geneID]

            pval = float(row[args.pval])

            logfc = float(row[args.logfc])

            if direction == "UP":
                if logfc < 0:
                    pval = 1

            elif direction == "DOWN":
                if logfc > 0:
                    pval = 1
            elif direction == "ANY":
                pass

            if pval < 0.05 and abs(logfc) > 1:
                significantGenes.add(geneID)
            else:
                measuredGenes.add(geneID)

        testSets = {}
        for setElem in availSets:

            setData = availSets[setElem]
            origGenes = setData[1]
            setGenes = [x for x in setData[1] if x in significantGenes or x in measuredGenes]

            removedGenes = origGenes.difference(setGenes)

            if len(removedGenes) > 0:
                remgen = ""

                if len(removedGenes) < 10:
                    remgen = ";".join(removedGenes)

                print(setElem, "removed genes: ", len(removedGenes), " ", remgen)

            testSets[setElem] = (setData[0], setGenes)

        from scipy.stats import hypergeom

        populationSize = len(significantGenes) + len(measuredGenes)
        numSuccInPopulation = len(significantGenes)

        setToResult = {}

        for setElem in testSets:
            setData = testSets[setElem]
            sampleSize = len(setData[1])
            successIntersection = significantGenes.intersection(setData[1])
            drawnSuccesses = len(successIntersection)

            pval = hypergeom.sf(drawnSuccesses - 1, populationSize, numSuccInPopulation, sampleSize)
            fractionOfHitSamples = drawnSuccesses / sampleSize if sampleSize > 0 else 0

            resultObj = {
                'elem_id': setElem,
                'population_size': populationSize,
                'success_population': numSuccInPopulation,
                'sample_size': sampleSize,
                'success_samples': drawnSuccesses,
                'pval': pval,
                'sample_success_fraction': fractionOfHitSamples,
                'genes': ";".join(successIntersection),
                'direction': direction
            }

            setToResult[setElem] = resultObj

        sortedElems = [x for x in setToResult]
        elemPvals = [setToResult[x]["pval"] for x in sortedElems]

        rej, elemAdjPvals, _, _ = multipletests(elemPvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

        for eidx, elem in enumerate(sortedElems):
            assert (setToResult[elem]['pval'] == elemPvals[eidx])
            setToResult[elem]['adj_pval'] = elemAdjPvals[eidx]

        for elem in sortedElems:
            dr = DataRow.fromDict(setToResult[elem])
            outdf.addRow(dr)

    outdf.export(args.output.name)
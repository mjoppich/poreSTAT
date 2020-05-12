

import argparse
import sys
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from collections import Counter

def loadEnhancement( fileE ):

    if fileE == None:
        print("Not loading gene name enhancements")
        return {}

    print("Loading gene name enhancements", fileE.name)

    ens2sym = {}

    for line in fileE:
        line = line.strip().split("\t")

        if not line[0].startswith("ENS"):
            continue

        ensemblID = line[0]
        geneSymbol = line[1]
        biotype = line[2]

        if len(geneSymbol) == 0:
            continue

        ens2sym[ensemblID] = (geneSymbol, biotype)

    return ens2sym

def loadGeneLengths(fileE):

    if fileE == None:
        print("Not loading gene lengths")
        return None

    print("Loading gene lengths", fileE.name)

    """
        Ensembl_gene_identifier GeneID  length
        ENSMUSG00000000001      14679   3262
        ENSMUSG00000000003      54192   902
        ENSMUSG00000000028      12544   2252
    """

    ens2gl = {}
    for lidx, line in enumerate(fileE):
        line = line.strip().split("\t")

        if lidx == 0:
            try:
                int(line[1])
            except:
                continue

        ensemblID = line[0]
        geneLength = line[1]

        if len(ensemblID) == 0 or len(geneLength) == 0:
            continue

        geneLength = int(geneLength)
        ens2gl[ensemblID] = geneLength

    return ens2gl

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-f', '--fc', type=argparse.FileType('r'), required=True, help='fc files')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True, help='output files')
    
    parser.add_argument('-e', '--enhanced', type=argparse.FileType('r'), default=None)
    parser.add_argument('-l', '--lengths', type=argparse.FileType('r'), required=True, default=None)
    parser.add_argument('-rrna', '--no-rrna', dest='norrna', action='store_true', default=False)

    parser.add_argument('-fpkm', '--fpkm', dest='fpkm', action='store_true', default=False)
    parser.add_argument('-tpm', '--tpm', dest='tpm', action='store_true', default=False)
    args = parser.parse_args()




    enhancedData = loadEnhancement(args.enhanced)
    geneLengths = loadGeneLengths(args.lengths)

    if args.norrna and enhancedData == None:
        raise argparse.ArgumentParser().error("removal of rRNA requires --enhanced!")

    indf = DataFrame.parseFromFile(args.fc.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
    })

    allheaders = indf.getHeader()
    featureCountsColumns = [y for y in ["Geneid", "Chr", "Start", "End", "Strand", "Length"] if y in allheaders]

    
    sampleHeaders = [x for x in allheaders if not x in featureCountsColumns]

    outdf = DataFrame()

    sample2total = Counter()
    sample2ratio = Counter()

    for row in indf:

        for sample in sampleHeaders:

            sample2total[sample] += row[sample]
            geneID = row["Geneid"]
            geneCount = row[sample]
            geneLength = geneLengths.get(geneID, 0)

            if geneLength != 0:
                sample2ratio[sample] += geneCount/geneLength

    allRowUpdates = []
    for row  in indf:

        curGeneID = row["Geneid"]

        if enhancedData != None and args.norrna:
            enhancedGeneData = enhancedData.get(curGeneID, None)
            
            # remove rRNA if needed
            if enhancedGeneData != None and "rRNA" in enhancedGeneData[1]:
                continue

        rowDict = {"Geneid": curGeneID}

        for x in featureCountsColumns:
            rowDict[x] = row[x]

        if not curGeneID in geneLengths:
            print("Missing gene ID", curGeneID)

        geneLength = geneLengths.get(curGeneID, 0)

        for sample in sampleHeaders:

            rowDict[sample] = row[sample]

            if args.fpkm:

                #print(curGeneID, geneLength)

                fpkmValue = row[sample]/(sample2total[sample]*geneLength) * pow(10,9)
                rowDict[sample + ".FPKM"] = fpkmValue

            if args.tpm:

                tpmValue = row[sample] / (geneLength * sample2ratio[sample]) * pow(10,6)
                rowDict[sample + ".TPM"] = tpmValue

        allRowUpdates.append(rowDict)

    allCols = set()
    for x in allRowUpdates:
        for y in x:
            if not y in featureCountsColumns:
                allCols.add(y)

    outdf.addColumns(featureCountsColumns)
    outdf.addColumns(sorted(allCols), default=0, ignoreDuplicates=True)
    outdf.updateRowIndexed("Geneid", allRowUpdates, ignoreMissingCols=True, addIfNotFound=True)

    outdf.export(args.output.name, exType=ExportTYPE.TSV)
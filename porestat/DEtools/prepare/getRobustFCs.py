

import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl
import argparse
import math

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d', '--de', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-p', '--cutoff', type=float, help='alignment files', default=0.05)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True, help='de files')

    args = parser.parse_args()


    for fidx, defile in enumerate(args.de):
        indf = DataFrame.parseFromFile(defile.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })

        inHeaders = indf.getHeader()

        genesymname = None

        if "gene_symbol" in inHeaders:
            genesymname = "gene_symbol"
        elif "Geneid" in inHeaders:
            genesymname = "Geneid"
        else:
            genesymname = "id"

        #if not args.gene in inHeaders:
        #    print("Unknown gene id column", args.gene)
        #    print(inHeaders)
        #    exit(-1)

        if not all(["ROB_log2FC" in inHeaders, "ROB_ADJ.PVAL" in inHeaders]): #"ROB_RAW.PVAL" in inHeaders,
            print("Not all necessary columns found.")
            exit(-1)



        for row in indf:

            if row["ROB_ADJ.PVAL"] < args.cutoff:

                rowVec = []
                for x in [genesymname, "ROB_log2FC", "ROB_RAW.PVAL", "ROB_ADJ.PVAL"]:
                    if x in row:
                        rowVec.append(row[x])
               
                print(*rowVec, sep="\t", file=args.output)


import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import matplotlib as mpl
import argparse
import math

import sys, os
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
import venn

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-p', '--pathways', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-n', '--top_n', type=int, required=False, default=100)
    parser.add_argument('-o', '--output', type=str, required=False, help="output base")


    args = parser.parse_args()


    if len(args.pathways) <= 2 or len(args.pathways) > 5:
        exit()



    pwcount2function = {
        2: venn.venn2,
        3: venn.venn3,
        4: venn.venn4,
        5: venn.venn5,
        6: venn.venn6,
    }


    foundRes = []

    for pathwaysFile in args.pathways:
        indf = DataFrame.parseFromFile(pathwaysFile.name, skipChar='#', replacements={
            "None": None,
            "": None,
            "NA": None
        })



        # try to find ID column

        idColumn = None
        inHeaders = indf.getHeader()

        for header in inHeaders:

            if "ID" in header:
                idColumn = header


        if idColumn == None:
            continue

        topNIDs = []
        for ridx, row in enumerate(indf):

            if ridx >= args.top_n:
                break



            topNIDs.append(row[idColumn])


        foundRes.append((pathwaysFile.name, topNIDs))


    vennLabels = venn.generate_petal_labels([set(x[1]) for x in foundRes])
    fig, ax = pwcount2function[len(foundRes)](vennLabels, names=[x[0] for x in foundRes])


    outname = args.output

    if not outname.endswith(".png"):
        outname += ".png"

    plt.savefig(outname, bbox_inches ="tight")


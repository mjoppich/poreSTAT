    

import argparse
import sys
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from collections import Counter
    
    
if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-f1', '--fc1', type=argparse.FileType('r'), required=True, help='fc files')
    parser.add_argument('-f2', '--fc2', type=argparse.FileType('r'), required=True, help='fc files')
    
    args = parser.parse_args()

    indf1 = DataFrame.parseFromFile(args.fc1.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
    })

    indf2 = DataFrame.parseFromFile(args.fc2.name, skipChar='#', replacements = {
        "None": None,
        "": None,
        "NA": None
    })


    featureCountsColumns = ["Geneid",	"Chr",	"Start",	"End",	"Strand",	"Length"]

    df1Samples = [x for x in indf1.getHeader() if not x in featureCountsColumns]
    df2Samples = [x for x in indf2.getHeader() if not x in featureCountsColumns]

    outdf = DataFrame()

    geneid2data = {}

    for row in indf1:

        data = {}

        for x in featureCountsColumns + df1Samples:

            data[x] = row[x]

        geneid2data[row["Geneid"]] = data


    for row in indf2:

        data = geneid2data.get(row["Geneid"], {})

        for x in df2Samples:
            data[x] = row[x]

        geneid2data[row["Geneid"]] = data


    outdf.addColumns(featureCountsColumns + df1Samples + df2Samples)

    allRowUpdates = [geneid2data[x] for x in sorted([y for y in geneid2data])]

    outdf.updateRowIndexed("Geneid", allRowUpdates, ignoreMissingCols=True, addIfNotFound=True)

    outdf.export("test1.tsv", ExportTYPE.TSV)
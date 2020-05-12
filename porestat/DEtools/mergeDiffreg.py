    

import argparse
import sys
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from collections import Counter
    
    
if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d1', '--de1', type=argparse.FileType('r'), required=True, help='de files')
    parser.add_argument('-d2', '--de2', type=argparse.FileType('r'), required=True, help='de files')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True, help='de files')

    parser.add_argument("-p1", "--prefix1", type=str, required=True)
    parser.add_argument("-p2", "--prefix2", type=str, required=True)

    parser.add_argument('-s2', '--samples', nargs='+', type=str, default=[])


    parser.add_argument('-pc', '--prefix-counts', dest="prefix_counts", action='store_true', default=False, help="run FC part")


    args = parser.parse_args()

    indf1 = DataFrame.parseFromFile(args.de1.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
    })

    indf2 = DataFrame.parseFromFile(args.de2.name, skipChar='#', replacements = {
        "None": None,
        "": None,
        "NA": None
    })

    allSamples = args.samples
    print("all samples", allSamples)

    df1Cols = indf1.getHeader()
    df2Cols = indf2.getHeader()

    df1UniqueCols = []
    df2UniqueCols = []

    df1SampleCols = []
    df2SampleCols = []
    
    for x in df1Cols:
        

        if x in args.samples:
            df1SampleCols.append(x)
            continue

        if x.startswith(tuple(args.samples)):
            df1SampleCols.append(x)
            continue

        if x in df2Cols:
            continue

        df1UniqueCols.append(x)

    for x in df2Cols:

        if x in args.samples:
            df2SampleCols.append(x)
            continue

        if x.startswith(tuple(args.samples)):
            df2SampleCols.append(x)
            continue

        if x in df1Cols:
            continue

        df2UniqueCols.append(x)

    df1SpecialCols = [x for x in df1Cols if any(["PVAL" in x, "FC" in x])]
    df2SpecialCols = [x for x in df2Cols if any(["PVAL" in x, "FC" in x])]

    df1UniqueCols = [x for x in df1UniqueCols if not x in df1SpecialCols]
    df2UniqueCols = [x for x in df2UniqueCols if not x in df2SpecialCols]

    df12CommonCols = [x for x in df1Cols if x in df2Cols and not x in df1SpecialCols and not x in df1SampleCols and not x in df2SampleCols and not x in df2SpecialCols]

    #df1SampleCols = [x for x in df1SampleCols if not x in df12CommonCols]
    #df2SampleCols = [x for x in df2SampleCols if not x in df12CommonCols]

    # any common column should contain the same information ...
    for ridx, (df1row, df2row) in enumerate(zip(indf1, indf2)):

        if ridx == 10:
            break


        allCommonCols = [x for x in df12CommonCols]
        for comCol in allCommonCols:

            if df1row[comCol] != df2row[comCol]:
                print("CHG", comCol)
                df12CommonCols.remove(comCol)
                df1UniqueCols.append(comCol)
                df2UniqueCols.append(comCol)


    df1UniqueCols = sorted(df1UniqueCols)
    df2UniqueCols = sorted(df2UniqueCols)
    df1SpecialCols = sorted(df1SpecialCols)
    df2SpecialCols = sorted(df2SpecialCols)
    df1SampleCols = sorted(df1SampleCols)
    df2SampleCols = sorted(df2SampleCols)
    df12CommonCols = sorted(df12CommonCols)

    df1Col2New = {}
    df2Col2New = {}


    for x in df1SampleCols:
        print("Sa1", x)
    for x in df2SampleCols:
        print("Sa2", x)

    for x in df12CommonCols:
        print("C", x)

    for x in df1UniqueCols:
        print("DF1", x)

    for x in df2UniqueCols:
        print("DF2", x)

    for x in df1SpecialCols:
        xn = x.split("_")
        xn.insert(1, args.prefix1)
        xn = "_".join(xn)
        df1Col2New[x] = xn
        print("Sp1", x, xn)
    for x in df2SpecialCols:
        xn = x.split("_")
        xn.insert(1, args.prefix2)
        xn = "_".join(xn)
        df2Col2New[x] = xn
        print("Sp2", x, xn)
    
    df1NewCols = [df1Col2New[x] for x in df1Col2New]
    df2NewCols = [df2Col2New[x] for x in df2Col2New]
    
    
    outdf = DataFrame()

    if args.prefix_counts:
        outdf.addColumns(df12CommonCols + [args.prefix1 + "_" + x for x in df1UniqueCols + df1SampleCols] + [args.prefix2 + "_" + x for x in df2UniqueCols+df2SampleCols] + df1NewCols + df2NewCols)
    else:
        outdf.addColumns(df12CommonCols + df1UniqueCols + df2UniqueCols + df1NewCols + df2NewCols)




    for x in outdf.getHeader():
        print("O", x)

    id2dataDf = {}

    for row in indf1:
        data = {}
        for x in df12CommonCols:
            data[x] = row[x]

        for x in df1UniqueCols+df1SampleCols:

            if args.prefix_counts:
                data[args.prefix1 + "_" +x] = row[x]
            else:
                data[x] = row[x]

        for x in df1Col2New:
            data[df1Col2New[x]] = row[x]

        id2dataDf[data["id"]] = data

    for row in indf2:
        data = id2dataDf.get(row["id"], {})

        for x in df12CommonCols:
            data[x] = row[x]

        for x in df2UniqueCols+df2SampleCols:
            if args.prefix_counts:
                data[args.prefix2 + "_" +x] = row[x]
            else:
                data[x] = row[x]


        for x in df2Col2New:
            data[df2Col2New[x]] = row[x]


    allRowUpdates = [id2dataDf[x] for x in sorted([y for y in id2dataDf])]

    outdf.updateRowIndexed("id", allRowUpdates, ignoreMissingCols=True, addIfNotFound=True)

    outdf.export(args.output.name, ExportTYPE.TSV)
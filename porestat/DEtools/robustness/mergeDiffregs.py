
import argparse
import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from collections import Counter


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-d', '--de', type=argparse.FileType('r'), nargs="+", required=True, help='de files')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True, help='de files')
    parser.add_argument("-p", "--prefixes", nargs='+', type=str, required=True)
    parser.add_argument('-s2', '--samples', nargs='+', type=str, default=[])

    parser.add_argument('-pc', '--prefix-counts', dest="prefix_counts", action='store_true', default=False, help="run FC part")


    args = parser.parse_args()

    curDF = DataFrame.parseFromFile(args.de[0].name, skipChar='#', replacements = {
        "None": None,
        "": None,
        "NA": None
    })

    for didx, deTable in enumerate(args.de):

        if didx == 0:
            continue

        indf2 = DataFrame.parseFromFile(deTable.name, skipChar='#', replacements = {
        "None": None,
        "": None,
        "NA": None
        })

        if didx == 1:
            curPrefix = args.prefixes[0]
        else:
            curPrefix = ""

        curCols = curDF.getHeader()
        df2Cols = indf2.getHeader()

        curUniqueCols = [x for x in curCols if (not x in df2Cols or x in args.samples)]
        df2UniqueCols = [x for x in df2Cols if (not x in curCols or x in args.samples)]

        curSpecialCols = [x for x in curCols if any(["PVAL" in x, "FC" in x])]
        df2SpecialCols = [x for x in df2Cols if any(["PVAL" in x, "FC" in x])]

        curUniqueCols = [x for x in curUniqueCols if not x in curSpecialCols]
        df2UniqueCols = [x for x in df2UniqueCols if not x in df2SpecialCols]

        curDF2CommonCols = [x for x in curCols if x in df2Cols and not x in curSpecialCols and not x in df2SpecialCols]

        # any common column should contain the same information ...
        for ridx, (df1row, df2row) in enumerate(zip(curDF, indf2)):

            if ridx == 10:
                break

            allCommonCols = [x for x in curDF2CommonCols]
            for comCol in allCommonCols:

                if df1row[comCol] != df2row[comCol]:
                    print("CHG", comCol)
                    curDF2CommonCols.remove(comCol)
                    curUniqueCols.append(comCol)
                    df2UniqueCols.append(comCol)

        df1UniqueCols = sorted(curUniqueCols)
        df2UniqueCols = sorted(df2UniqueCols)
        df1SpecialCols = sorted(curSpecialCols)
        df2SpecialCols = sorted(df2SpecialCols)
        df12CommonCols = sorted(curDF2CommonCols)

        df1Col2New = {}
        df2Col2New = {}

        for x in curDF2CommonCols:
            print("C", x)
        for x in curUniqueCols:
            print("DF1", x)
        for x in df2UniqueCols:
            print("DF2", x)
        for x in curSpecialCols:
            xn = x.split("_")

            if len(curPrefix) > 0:
                xn.insert(1, curPrefix)

            xn = "_".join(xn)
            df1Col2New[x] = xn
            print("S1", x, xn)
        for x in df2SpecialCols:
            xn = x.split("_")
            xn.insert(1, args.prefixes[didx])
            xn = "_".join(xn)
            df2Col2New[x] = xn
            print("S2", x, xn)

        df1NewCols = [df1Col2New[x] for x in df1Col2New]
        df2NewCols = [df2Col2New[x] for x in df2Col2New]

        outdf = DataFrame()

        if args.prefix_counts:

            if len(curPrefix) > 0:
                curPrefix += "_"

            outdf.addColumns(df12CommonCols + [curPrefix + x for x in df1UniqueCols] + [args.prefixes[didx] + "_" + x for x in df2UniqueCols] + df1NewCols + df2NewCols)
        else:
            outdf.addColumns(df12CommonCols + df1UniqueCols + df2UniqueCols + df1NewCols + df2NewCols)

        for x in outdf.getHeader():
            print("O", x)

        id2dataDf = {}

        for row in curDF:
            data = {}
            for x in df12CommonCols:
                data[x] = row[x]

            for x in df1UniqueCols:

                if args.prefix_counts:
                    data[args.prefix1 + "_" + x] = row[x]
                else:
                    data[x] = row[x]

            for x in df1Col2New:
                data[df1Col2New[x]] = row[x]

            id2dataDf[data["id"]] = data

        for row in indf2:
            data = id2dataDf.get(row["id"], {})

            for x in df12CommonCols:
                data[x] = row[x]

            for x in df2UniqueCols:
                if args.prefix_counts:
                    data[args.prefix2 + "_" + x] = row[x]
                else:
                    data[x] = row[x]

            for x in df2Col2New:
                data[df2Col2New[x]] = row[x]

        allRowUpdates = [id2dataDf[x] for x in sorted([y for y in id2dataDf])]

        outdf.updateRowIndexed("id", allRowUpdates, ignoreMissingCols=True, addIfNotFound=True)
        curDF = outdf

    curDF.export(args.output.name, ExportTYPE.TSV)
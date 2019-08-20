import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--detable', nargs='+', type=argparse.FileType('r'), help='de table file to read in',
                        required=True)
    parser.add_argument('-i', '--include', type=str, nargs='+', help="Location to store biogui templates in", default=[],
                        required=False)
    parser.add_argument('-e', '--exclude', type=str, nargs='+', help="prefix for all variables", default=[], required=False)
    parser.add_argument('-c', '--conf-level', type=float, help="confidence level for pvalue", default=0.05)
    parser.add_argument('-v', '--verbose', help="verbose output", required=False, default=False, action='store_true')
    parser.add_argument('-fc', '--fc', help="gene,foldchange", required=False, default=False, action='store_true')

    args = parser.parse_args()

    for detable in args.detable:

        containedMethods = []
        method2colidx = {}
        method2fc = {}

        for idx, line in enumerate(detable):

            aline = line.strip().split("\t")

            if idx == 0:

                # DESeq_log2FC	DESeq_RAW.PVAL	DESeq_ADJ.PVAL	msEmpiRe_log2FC	msEmpiRe_RAW.PVAL	msEmpiRe_ADJ.PVAL

                for elem in aline:
                    if elem.endswith("ADJ.PVAL"):
                        aelem = elem.split("_")
                        aelem = "_".join(aelem[0:len(aelem) - 1])

                        containedMethods.append(aelem)

                #print("Found DE methods:", containedMethods)

                for cidx, elem in enumerate(aline):

                    for method in containedMethods:

                        if elem == method + "_ADJ.PVAL":
                            method2colidx[method] = cidx

                        if elem == method + "_log2FC":
                            method2fc[method] = cidx

                """
                for x in containedMethods:

                    if not x in method2colidx:
                        print("No col idx found for method:", x)

                for x in method2colidx:
                    print(x, method2colidx[x])
                    
                """

                for x in args.include + args.exclude:
                    assert(x in method2colidx)

            else:

                allArgsLessConf = True
                for x in args.exclude:

                    if x in method2colidx:

                        melem = aline[method2colidx[x]]

                        if not (melem == 'None' or melem == "NA") and float(melem) < args.conf_level:
                            allArgsLessConf = False

                if not allArgsLessConf:
                    continue

                allArgsConf = True
                includeFCs = set()
                for x in args.include:

                    if x in method2colidx:

                        melem = aline[method2colidx[x]]

                        if (melem == 'None' or melem == "NA") or float(melem) >= args.conf_level:
                            allArgsConf = False
                        else:
                            includeFCs.add(float(aline[method2fc[x]]))

                if not allArgsConf:
                    continue

                if args.fc:
                    print(aline[0], min(includeFCs), sep="\t")
                else:
                    print(aline[0], sep="\t")

                if args.verbose:
                    for x in method2colidx:
                        print(x, aline[method2colidx[x]], aline[method2fc[x]])
import argparse
from collections import defaultdict

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list', nargs='+', type=argparse.FileType('r'), help='de table file to read in',
                        required=True)
    parser.add_argument('-m', '--mode', type=str, help="prefix for all variables", default="MIN", required=False)
    parser.add_argument('-s', '--select', type=str, help="prefix for all variables", default="ALL", required=False)

    args = parser.parse_args()

    curMode = args.mode.upper()
    curSelect = args.select.upper()

    if not curMode in ["MIN", "AVG"]:
        parser.error("MODE MUST BE IN MIN,AVG")

    if not curSelect in ["UNIQ", "INTER", "ALL"]:
        parser.error("SELECT MUST BE IN UNIQ, INTER, ALL")

    gene2fc = defaultdict(list)

    for infile in args.list:

        for line in infile:

            line = line.strip().split("\t")

            gene = line[0]
            fc = float(line[1])

            gene2fc[gene].append(fc)


    for gene in sorted([x for x in gene2fc]):

        fcs = gene2fc[gene]

        fc = None
        if curMode == 'MIN':

            fc = min(fcs)

        else:
            fc = sum(fcs)/len(fcs)

        if fc == None:
            continue

        if curSelect == "ALL":
            print(gene, fc, sep="\t")
        elif curSelect == "UNIQ":
            if len(fcs) == 1:
                print(gene, fc, sep="\t")
        elif curSelect == "INTER":
            if len(fcs) == len(args.list):
                print(gene, fc, sep="\t")


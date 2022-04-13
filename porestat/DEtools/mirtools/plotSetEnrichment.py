from collections import defaultdict

import matplotlib.pyplot as plt

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")

import numpy as np
import argparse

def transformPVal(pval):
    return -np.log2(pval)


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-p', '--pvalue', type=int, default=7)
    parser.add_argument('-n', '--set-name', type=int, default=0)
    parser.add_argument('-t', '--tsv', type=argparse.FileType('r'), nargs="+", required=True, help='de files')
    args = parser.parse_args()



    for inputfile in args.tsv:

            """

        elem_id	population_size	success_population	sample_size	success_samples	sample_success_fraction	pval	adj_pval	genes
        histone	38058	2911	521	51	0.09788867562380038	0.04224174792964982	0.21120873964824907	ERBB4;GDNF;KDM6B;ARNTL;GET1;FOXO1;KDM6A;GATAD2A;SPOP;CBX7;PRDM5;PRDM6;PRKAA2;EZH2;JMJD1C;ATXN7;UBE2D1;DOT1L;SMYD1;HMGN3;ASXL3;KANSL2;CBX1;RNF20;TLE4;KDM5C;PRKAG3;USP44;PPP4C;VDR;DNAJC2;USP36;HMGN5;DNTTIP2;ASXL1;GADD45G;CBX5;BRD1;HCFC2;FOXP2;NAP1L2;ZBTB16;GADD45B;TAF9B;HDAC9;SRCAP;ZNF711;JARID2;BANP;CDK17;PADI4
        rna	38058	2911	33	1	0.030303030303030304	0.9277069476710338	1.0	DDX21
        chromatin	38058	2911	75	5	0.06666666666666667	0.688470584211667	1.0	ACTL6B;RNASE12;CHD7;VPS72;PSIP1
        dna	38058	2911	100	6	0.06	0.7857110464381157	1.0	GDNF;CHD1;FOXO1;FOXP2;CHD2;TET1
        tf	38058	2911	1	0	0.0	1.0	1.0	

            """

            import numpy as np

            foundSets = defaultdict(list)

            hasDirection = False

            allElemNames = set()
            for lidx, line in enumerate(inputfile):

                if lidx == 0:

                    if "direction" in line:
                        hasDirection = True

                    continue

                line = line.split("\t")

                elemName = line[args.set_name]
                allElemNames.add(elemName)

                pval = float(line[ args.pvalue ])

                if hasDirection:
                    foundSets[elemName].append((pval, line[8], int(line[4])))


            elem2sigdir = {}

            for elemName in foundSets:

                elemRes = foundSets[elemName]

                if len(elemRes) < 2:
                    seenDirs = set([x[1] for x in elemRes])

                    if not "UP" in seenDirs:
                        elemRes.append((1.0, "UP", 1))

                    if not "DOWN" in seenDirs:
                        elemRes.append((1.0, "DOWN", 1))


                if len(elemRes) > 1:

                    countSig = sum([1 for x in elemRes if x[0] < 0.05])
                    SIGDIR = [x[1] for x in elemRes if x[0] < 0.05]
                    elem2sigdir[elemName] = (SIGDIR, tuple(sorted(elemRes, key=lambda x: x[1])))

                        # print(elemName, SIGDIR,  "; ".join([str(x) for x in elemRes]))


            def getElemScore(elem):

                logUp = transformPVal(elem[-1][0])
                logDown = transformPVal(elem[-2][0])

                maxScore = max(logUp, logDown)

                logAny = 0
                if len(elem) > 2:
                    logAny = transformPVal(elem[0][0])

                pValDiff = abs(logUp-logDown)

                if logUp < transformPVal(0.05) and logDown < transformPVal(0.05):
                    pValDiff = 0
                    maxScore = -1000

                if logUp > transformPVal(0.05) and logDown > transformPVal(0.05):
                    pValDiff = 0
                    maxScore = -1000

                anyDist = abs(max(logUp, logDown)-logAny)

                return maxScore, pValDiff


            elemKeys = sorted([x for x in elem2sigdir], reverse=True, key=lambda x: getElemScore(elem2sigdir[x][1]))

            myMirs = ["miR_29a", "miR_29b", "miR_29c"]
            accmyMirs = []
            for x in myMirs:
                print(x, elemKeys.index(x) if x in elemKeys else -1)
                if x in elemKeys:
                    accmyMirs.append(x)


            maxElems = int(sys.argv[2]) if len(sys.argv) >= 3 else 10

            if len(elemKeys) > maxElems:
                foundElems = elemKeys[0:maxElems]
            else:
                foundElems = elemKeys



            elemLabels = []
            elemValues = []

            for foundElem in foundElems+accmyMirs:

                elemData = elem2sigdir[foundElem]

                dir2perf = {}
                for x in elemData[1]:
                    dir2perf[x[1]] = x[0]

                for dir in ["ANY", "DOWN", "UP"]:
                    elemLabels.append(foundElem.replace("_", "-") + " {}".format(dir))
                    elemValues.append(dir2perf.get(dir, 1))

            elemValues = [transformPVal(x) for x in elemValues]

            print(elemValues[:5])
            print(elemLabels[:5])

            if len(elemValues) > 1:

                y_pos = np.arange(len(elemLabels))
                fig, ax = plt.subplots(figsize=(8, 5+maxElems*0.5))

                ax.barh(y_pos, elemValues, align='center')
                ax.set_yticks(y_pos)
                ax.set_yticklabels(elemLabels)
                ax.invert_yaxis()  # labels read top-to-bottom
                ax.set_xlabel('-log2(adj_pval)')
                ax.set_title('Top enriched Gene-Sets')

                ax.axvline(x=transformPVal(0.05), ymin=y_pos[0], ymax=y_pos[-1], c="r")

                outname = inputfile.name + ".png"
                plt.savefig(outname, bbox_inches="tight")
                print(outname)
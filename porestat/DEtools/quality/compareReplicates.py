

import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl
import argparse
import matplotlib.patches as patches
import sys, os
import math
sys.path.insert(0, "/mnt/d/dev/git/poreSTAT/")

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from collections import Counter
import natsort

mpl.style.use("seaborn")




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-d', '--conditions', nargs='+', type=str, action='append',  required=True, help='alignment files')
    parser.add_argument('-p', '--pathname', action="store_true", default=False)
    parser.add_argument('-o', '--output', nargs='+', type=str, required=False, help="output base")

    parser.add_argument('-r', '--relative', action="store_true", default=False)
    args = parser.parse_args()

    if args.output == None:
        args.output = [counts.name for counts in args.counts]




    for fidx, defile in enumerate(args.counts):
        indf = DataFrame.parseFromFile(defile.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })

        inHeaders = indf.getHeader()

        #if not args.gene in inHeaders:
        #    print("Unknown gene id column", args.gene)
        #    print(inHeaders)
        #    exit(-1)

        for conditions in args.conditions:

            for condition in conditions:

                if not condition in inHeaders:
                    print("Unknown condition", condition)
                    print("Known conditions", inHeaders)
                    exit(-1)

        for cidx, conditions in enumerate(args.conditions):

            fig, axes = plt.subplots(len(conditions), len(conditions), figsize=(24, 18), sharex=True, sharey=True)
            fig.suptitle("Comparison of read counts", fontsize=15)

            if conditions[0].startswith("X."):
                sconditions = [".".join(x.split(".")[-4:-2]) for x in conditions]
            else:

                if args.pathname:
                    # just in case it's not an actual path :D
                    sconditions = []
                    for x in conditions:
                        nx = x.split(os.path.sep)
                        sconditions.append(nx[max([0, len(nx)-2])])

                else:
                    sconditions = [os.path.basename(x) for x in conditions]

            sconditions = ["\n".join(x.split(".")) for x in sconditions]
                
            axes[0,0].set_title(sconditions[0])

            if args.relative:
                print("condition1", "condition2", "condition1_total", "condition2_total", "scaling_factor_cond1", "scaling_factor_cond2", "scaling_factor", "union_gene_count", "intersect_gene_count", "intersect_gene_count_fraction1", "intersect_gene_count_fraction2", "least_frequent_cond1", "least_frequent_cond2")

            cond2genecount = {}

            for i in range(0, len(conditions)):

                # build a rectangle in axes coords
                left, width = .25, .5
                bottom, height = .25, .5
                right = left + width
                top = bottom + height

                # axes coordinates are 0,0 is bottom left and 1,1 is upper right
                p = patches.Rectangle(
                    (left, bottom), width, height,
                    fill=False, transform=axes[i,i].transAxes, clip_on=False
                    )

                axes[i,i].add_patch(p)

                axes[i,i].text(0.5*(left+right), 0.5*(bottom+top), sconditions[i],
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=7, color='red',
                transform=axes[i,i].transAxes)

                for j in range(i+1, len(conditions)):

                    cond1counts = Counter()
                    cond2counts = Counter()

                    cond1countsRaw = Counter()
                    cond2countsRaw = Counter()

                    cond1totals = 0
                    cond2totals = 0

                    for row in indf:

                        cond1count = row[conditions[i]]
                        cond2count = row[conditions[j]]

                        gene = row["Geneid"]

                        cond1countsRaw[gene] = cond1count
                        cond2countsRaw[gene] = cond2count

                        # required for log!
                        if cond1count == 0:
                            cond1count = 0.01

                        # required for log!
                        if cond2count == 0:
                            cond2count = 0.01

                        cond1totals += cond1count
                        cond2totals += cond2count

                        

                        cond1counts[gene] = cond1count
                        cond2counts[gene] = cond2count

                    cond2genecount[conditions[i]] = cond1countsRaw
                    cond2genecount[conditions[j]] = cond2countsRaw

                    plotGenes = set([x for x in cond1counts]).intersection([x for x in cond2counts])

                    if args.relative:
                        sf1 = math.pow(10, math.ceil(math.log10(cond1totals)))
                        sf2 = math.pow(10, math.ceil(math.log10(cond2totals)))
                        plotGenes = set([x[0] for x in cond1counts.most_common(2000)] + [x[0] for x in cond2counts.most_common(2000)])
                        plotGenesInt = set([x[0] for x in cond1counts.most_common(2000)]).intersection([x[0] for x in cond2counts.most_common(2000)])

                        cond1Top = sum([cond1counts[x] for x in plotGenesInt]) / cond1totals
                        cond2Top = sum([cond2counts[x] for x in plotGenesInt]) / cond2totals

                        cond1Least = "__".join([str(x) for x in cond1counts.most_common(2000)[-1]])
                        cond2Least = "__".join([str(x) for x in cond2counts.most_common(2000)[-1]])
                        
                        sf = max([sf1, sf2])
                        sf = 10000

                        print(conditions[i], conditions[j], cond1totals, cond2totals, sf1, sf2, sf, len(plotGenes), len(plotGenesInt), cond1Top, cond2Top, cond1Least, cond2Least)

                    else:
                        sf = 1.0
                        cond1totals = 1.0
                        cond2totals = 1.0



                    cond1countsNormed = [sf * cond1counts[x] / cond1totals for x in plotGenes]
                    cond2countsNormed = [sf * cond2counts[x] / cond2totals for x in plotGenes]

                    axes[i,j].scatter(cond1countsNormed, cond2countsNormed, label="{} vs. {} (n={})".format(sconditions[i], sconditions[j], len(plotGenes)) ,s=5)
                    axes[i,j].set_yscale('log')
                    axes[i,j].set_xscale('log')

                
                    if i == 0:
                        axes[i,j].set_title(sconditions[j])
                
            plt.savefig(args.output[fidx] + ".replicates."+str(fidx) + "." + str(cidx) +".png", bbox_inches ="tight")
            plt.close()

            topGeneCounts = [2000, 3000, 4000, 5000]
            headers = ["Samples"]
            for x in topGeneCounts:
                headers += ["Top{}%".format(x), "Count@{}".format(x)]
            headers += ["Genes >= 10", "Gene at border 10"]
            headers += ["Genes >= 1", "Gene at border 1"]
            headers += ["Total Counts"]

            print("\t".join(headers))
            

            for cond in natsort.natsorted(cond2genecount):

                cnter = cond2genecount[cond]

                total = sum([cnter[x] for x in cnter])

                repPercs = []
                for tg in topGeneCounts:

                    mcgenes = cnter.most_common(tg)
                    mccount = sum([x[1] for x in mcgenes])
                    lccount = mcgenes[-1][1]

                    repPercs.append(str(mccount/total))
                    repPercs.append(str(lccount))

                cnt10Idx = 0
                cnt10Name = None
                for gn, gc in cnter.most_common():
                    if gc < 10:
                        break

                    cnt10Idx += 1
                    cnt10Name = (gn, gc)

                repPercs.append(str(cnt10Idx))
                repPercs.append(str(cnt10Name))

                cnt0Idx = 0
                cnt0Name = None
                for gn, gc in cnter.most_common():
                    if gc < 1:
                        break

                    cnt0Idx += 1
                    cnt0Name = (gn, gc)

                repPercs.append(str(cnt0Idx))
                repPercs.append(str(cnt0Name))

                totalCount = sum([cnter[x] for x in cnter])
                repPercs.append(str(totalCountgit))

                print(cond, "\t".join(repPercs), sep="\t")
                    
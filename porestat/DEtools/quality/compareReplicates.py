

import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl
import argparse
import matplotlib.patches as patches
import math

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from collections import Counter
import natsort
from scipy.stats import pearsonr

mpl.style.use("seaborn-v0_8")




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-d', '--conditions', nargs='+', type=str, action='append',  required=True, help='alignment files')
    parser.add_argument('-p', '--pathname', action="store_true", default=False)
    parser.add_argument('-o', '--output', nargs='+', type=str, required=False, help="output base")

    parser.add_argument('-ls', '--ls', action="store_true", default=False )

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

        for conditions in args.conditions:

            for condition in conditions:

                if not condition in inHeaders:
                    print("Unknown condition", condition)
                    print("Known conditions", inHeaders)
                    exit(-1)

        for cidx, rawConditions in enumerate(args.conditions):

            batchSize = 20

            print("Creating Batches for elements:", len(rawConditions))

            conditionBatches = []
            if len(rawConditions) < batchSize:
                conditionBatches.append(rawConditions)
            else:
                curStart = 0
                while curStart < len(rawConditions):
                    conditionBatches.append(rawConditions[curStart:curStart+batchSize])
                    print(curStart, curStart+batchSize)
                    curStart += batchSize
            print("Got Batches:", len(conditionBatches))
            
            for batchIdx, conditions in enumerate(conditionBatches):
                
                print("Starting Batch", batchIdx)
                numConditions = len(conditions)

                maxCount = 1
                minCount = 1
                
                for condition in conditions:
                    condCount = sum(indf.getColumn(condition))
                    for row in indf:
                        if not args.ls:
                            maxCount = max(row[condition]+1, maxCount)
                            minCount = min(row[condition]+1, minCount)
                        else:
                            curCount = row[condition]/condCount
                            maxCount = max(curCount, maxCount)

                            if curCount > 0:
                                minCount = min(curCount, minCount)

                print("Min Count", minCount)
                print("Max Count", maxCount)

                print("Min Count", minCount, math.floor(math.log10(minCount)))
                print("Max Count", maxCount, math.ceil(math.log10(maxCount)))

                plotRange = (10** math.floor(math.log10(minCount)), 10**math.ceil(math.log10(maxCount)))
                print("Log Range", plotRange )
                plotRangeLog = (math.log10(plotRange[0]), math.log10(plotRange[1]))

                print("Starting figure for conditions", numConditions)

                fig, axes = plt.subplots(numConditions, numConditions, figsize=(max(20, numConditions*4), max(20, numConditions*4)))#, sharex=True, sharey=True)
                figTitle = "Comparison of Read Counts"
                if args.ls:
                    figTitle += " (library-size normalized)"
                fig.suptitle(figTitle, fontsize=24)

                print("Created Subplots")

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
                    
                

                cond2genecount = {}

                for i in range(0, len(conditions)):

                    # build a rectangle in axes coords
                    left, width = .25, .5
                    bottom, height = .25, .5
                    right = left + width
                    top = bottom + height

                    # axes coordinates are 0,0 is bottom left and 1,1 is upper right
                    #p = patches.Rectangle(
                    #    (left, bottom), width, height,
                    #    fill=False, transform=axes[i,i].transAxes, clip_on=False
                    #    )
                    #axes[i,i].add_patch(p)
                    #axes[i,i].text(0.5*(left+right), 0.5*(bottom+top), sconditions[i],
                    #horizontalalignment='center',
                    #verticalalignment='center',
                    #fontsize=7, color='red',
                    #transform=axes[i,i].transAxes)

                    axes[i,i].remove()
                    axes[i,i] = fig.add_subplot(len(conditions), len(conditions), (i*len(conditions)+i)+1)

                    if not args.ls:
                        cond1Counts = [x+1 for x in indf.getColumn(conditions[i])]
                    else:
                        condCount = sum(indf.getColumn(condition))
                        cond1Counts = [x/condCount for x in indf.getColumn(conditions[i])]

                    axes[i,i].hist(cond1Counts, bins=len(cond1Counts), histtype='step', density=True, cumulative=True,linewidth=2)
                    axes[i,i].set_xlim(plotRange)
                    axes[i,i].set_xscale("log")


                    for j in range(i+1, len(conditions)):

                        print(conditions[i], "vs.", conditions[j])

                        cond1counts = Counter()
                        cond2counts = Counter()

                        cond1countsRaw = Counter()
                        cond2countsRaw = Counter()

                        cond1totals = 0
                        cond2totals = 0

                        for row in indf:

                            cond1count = row[conditions[i]]
                            cond2count = row[conditions[j]]

                            gene = row[genesymname]

                            cond1countsRaw[gene] = cond1count
                            cond2countsRaw[gene] = cond2count

                            cond1totals += cond1count
                            cond2totals += cond2count

                            cond1counts[gene] = cond1count
                            cond2counts[gene] = cond2count

                        cond2genecount[conditions[i]] = cond1countsRaw
                        cond2genecount[conditions[j]] = cond2countsRaw

                        plotGenes = set([x for x in cond1counts]).intersection([x for x in cond2counts])

                        cond1countsIntersected = [cond1counts[x] for x in plotGenes]
                        cond2countsIntersected = [cond2counts[x] for x in plotGenes]

                        if args.ls:
                            cond1countsIntersected = [cond1counts[x]/cond1totals for x in plotGenes]
                            cond2countsIntersected = [cond2counts[x]/cond2totals for x in plotGenes]


                        if len(conditions) > 5:
                            print("More than 5 conditions => hexbin", len(conditions))

                            xMinPosValue = min([x for x in cond1countsIntersected if x > 0])
                            yMinPosValue = min([x for x in cond2countsIntersected if x > 0])

                            axes[i,j].hexbin([x+xMinPosValue for x in cond1countsIntersected], [x+yMinPosValue for x in cond2countsIntersected], label="{} vs. {} (n={})".format(sconditions[i], sconditions[j], len(plotGenes)), gridsize=25, cmap="plasma", xscale="log", yscale="log", bins="log", extent=[plotRangeLog[0], plotRangeLog[1],plotRangeLog[0], plotRangeLog[1]])

                        else:
                            print("LEQ than 5 conditions => hexbin", len(conditions))
                            axes[i,j].scatter(cond1countsIntersected, cond2countsIntersected, label="{} vs. {} (n={})".format(sconditions[i], sconditions[j], len(plotGenes)) ,s=5)
                            axes[i,j].set_yscale('log')
                            axes[i,j].set_xscale('log')




                        axes[i,j].plot([plotRange[0], plotRange[1]], [plotRange[0], plotRange[1]], color="red", linestyle="dotted")
                        axes[i,j].set_xlim(plotRange)
                        axes[i,j].set_ylim(plotRange)


                        # calculate correlation!
                        corrEff, corrSig = pearsonr(cond1countsIntersected, cond2countsIntersected)

                        corrText = "Corr: {:.3f}\nSig: {:.3f}".format(corrEff, corrSig)

                        p = patches.Rectangle(
                            (left, bottom), width, height,
                            fill=False, transform=axes[i,i].transAxes, clip_on=False
                            )
                        axes[j,i].add_patch(p)
                        axes[j,i].text(0.5*(left+right), 0.5*(bottom+top), corrText,
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=24, color='red',
                        transform=axes[j,i].transAxes)

                        axes[j,i].axes.xaxis.set_ticklabels([])
                        axes[j,i].axes.yaxis.set_ticklabels([])
                        axes[j,i].grid(False)

                    
                        if i == 0:
                            axes[i,j].set_title(sconditions[j])
                    

                axes[0,0].set_title(sconditions[0])

                plotFileName = args.output[fidx] + ".replicates."+str(fidx) + "." + str(cidx) + "." + str(batchIdx) +".png"
                print("Saving Plot", plotFileName)
                plt.savefig(plotFileName, bbox_inches ="tight")
                plt.close()
                print("Saving Plot Done", plotFileName)


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
                    repPercs.append(str(totalCount))

                    print(cond, "\t".join(repPercs), sep="\t")
                        






#python3 /mnt/f/dev/git/poreSTAT/porestat/DEtools/quality/compareReplicates.py --pathname --counts rawdata/counts/counts.reads.exon.all.counts.tsv --conditions k4_0_ccs_ccs k2_0_ccs_ccs k30_0_ccs_ccs k5_0_ccs_ccs k15_0_ccs_ccs k10_0_ccs_ccs k23_0_ccs_ccs --conditions m14_3_acs_acs_noinf m26_2_acs_acs_noinf m14_4_acs_acs_noinf m3_3_acs_acs_noinf m24_4_acs_acs_sub m25_2_acs_acs_noinf m21_1_acs_acs_inf --output diffregs/STEMI_PMN_CCS_ACS/STEMI_PMN_CCS_ACS.readsexon.diffreg/orig_counts.countreplicates
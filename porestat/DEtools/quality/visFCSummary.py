import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl
import argparse
import math

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../")


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


#mpl.style.use("seaborn")




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-s', '--summary', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-o', '--output', type=str, required=False, help="output base")

    args = parser.parse_args()

    if args.output == None:
        args.output = args.summary.name

    indf = DataFrame.parseFromFile(args.summary.name, skipChar='#', replacements = {
            "None": None,
            "": None,
            "NA": None
        })


    allStatus = []

    allCols = indf.getHeader()
    allCols.remove("Status")

    for row in indf:
        allStatus.append(row["Status"])

    sampleData = defaultdict(lambda: dict())

    for row in indf:
        for sample in allCols:       
            sampleData[sample][row["Status"]] = int(row[sample])


    for sample in sampleData:

        tc = 0
        for stat in sampleData[sample]:

            tc += sampleData[sample][stat]

        sampleData[sample]["Total"] = tc

    allStatus.insert(0, "Total")

    def autolabel(ax, rects):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                    str(round(height, 2)),
                    ha='center', va='bottom')



    for sidx, sample in enumerate(sampleData):

        ind = range(0, len(allStatus))
        
        fig, ax = plt.subplots()
        
        allValues = [sampleData[sample].get(x, 0) for x in allStatus]
        rects = ax.bar(ind, allValues)
        autolabel(ax, rects)

        plt.title(sample)

        plt.xticks(ind, allStatus,rotation=70)
        
        plt.xlabel("Category Count")
        plt.ylabel("Assigned Count")

        plt.savefig(args.output + ".fcsummary." + str(sidx) + ".png", bbox_inches ="tight")
        plt.close()


        ind = range(0, len(allStatus))
        fig, ax = plt.subplots()

        allValues = [sampleData[sample].get(x, 0)/sampleData[sample]["Total"] for x in allStatus]

        rects = ax.bar(ind, allValues)
        autolabel(ax, rects)

        plt.title(sample)

        plt.xticks(ind, allStatus,rotation=70)
        
        plt.xlabel("Category Rel Count")
        plt.ylabel("Assigned Rel Count")

        plt.savefig(args.output + ".fcsummary.rel." + str(sidx) + ".png", bbox_inches ="tight")
        plt.close()
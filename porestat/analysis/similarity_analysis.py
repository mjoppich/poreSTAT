import argparse
import HTSeq
import matplotlib

from porestat.plots.poreplot import PorePlot, MultiAxesPointHTMLTooltip

from porestat.plots.plotconfig import PlotConfig
from ..utils.DataFrame import DataFrame, DataRow
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Files import fileExists
from scipy import stats
import sys, math
import numpy as np
from matplotlib import pyplot as plt
import mpld3
import random
import os


class SimilarityAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(SimilarityAnalysisFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('similarity', help='expls help')
        parser.add_argument('-c', '--counts', nargs='+', type=str, help='counts summary file', required=False)

        def fileOpener( filename ):
            open(filename, 'w').close()
            return filename

        parser.add_argument('-o', '--output', type=fileOpener, help='output location, default: std out', default=sys.stdout)

        parser = PlotConfig.addParserArgs(parser)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return SimilarityAnalysis(simArgs)



class SimilarityAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(SimilarityAnalysis, self).__init__(args)
        self.counts = None

    def _makePropDict(self):

        return None

    def readCounts(self, args):

        for countsFile in args.counts:

            if not fileExists(countsFile):
                raise PSToolException("Read info file does not exist: " + str(countsFile))

        self.counts = {}

        for countsFile in args.counts:
            self.counts[countsFile] = DataFrame.parseFromFile(countsFile)

    def prepareInputs(self, args):
        self.readCounts(args)

        self.writeLinesToOutput(args.output, "\t".join(['gene', 'coverage', 'rank']) + "\n", mode='w')

        return []

    def execParallel(self, data, environment):

        return None


    def joinParallel(self, existResult, newResult, oEnvironment):

        return None


    def makeResults(self, parallelResult, oEnvironment, args):

        genes = {}
        unionGenes = set()

        for x in self.counts:
            df = self.counts[x]
            allgenes = df['gene']

            dfGenes = allgenes.to_set()

            unionGenes = unionGenes.union(dfGenes)

            genes[x] = dfGenes

        N = len(self.counts)
        unionGenes = list(unionGenes)

        dataRows = [x for x in self.counts]

        cntData = np.zeros((N,len(unionGenes)))

        for i in range(0, N):

            expName = dataRows[i]
            df = self.counts[expName]

            for j in range(0,len(unionGenes)):

                genename = unionGenes[j]
                geneExpr = df.findRow('gene', genename)

                if geneExpr == None:
                    cntData[i, j] = 0.0
                else:
                    cntData[i, j] = geneExpr['coverage']


        tooltipCSS = """
        table
        {
          border-collapse: collapse;
        }
        th
        {
          color: #ffffff;
          background-color: #000000;
        }
        td
        {
          background-color: #cccccc;
        }
        table, th, td
        {
          font-family:Arial, Helvetica, sans-serif;
          border: 1px solid black;
          text-align: right;
        }
        """

        def makeHTMLTable( gene, cond1, cond1val, cond2, cond2val):

            cond1 = cond1.split( os.path.sep )[-2]
            cond2 = cond2.split(os.path.sep)[-2]

            lines = []
            lines.append("<tr><th></th><th>{gn}</th></tr>".format(gn=gene))
            lines.append( "<tr><td>{c}</td><td>{cv:.2f}</td></tr>".format(c=cond1, cv=cond1val) )
            lines.append( "<tr><td>{c}</td><td>{cv:.2f}</td></tr>".format(c=cond2, cv=cond2val))



            htmlStr = "<table>"+"".join(lines)+"</table>"

            return htmlStr

        dataShape = cntData.shape

        fig, ax = plt.subplots(dataShape[0], dataShape[0], sharex="col", sharey="row", figsize=(8, 8))
        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
                            hspace=0.1, wspace=0.1)

        minx=float('inf')
        maxx=float('-inf')
        miny=float('inf')
        maxy=float('-inf')

        allPltAnnot = []
        for i in range(dataShape[0]):
            for j in range(dataShape[0]):
                xdata = []
                ydata = []
                colors = []

                X = cntData[j, :]
                Y = cntData[i, :]

                htmlData = []
                for k in range(0, len(unionGenes)):
                    #xdata.append(x[0])
                    #ydata.append(x[1])

                    if i != j and X[k] == Y[k] and X[k] > 0:
                        print(unionGenes[k])

                    htmlData.append( makeHTMLTable(unionGenes[k], dataRows[j], X[k], dataRows[i], Y[k]) )
                    colors.append( 'blue' )

                minx = min(minx, np.amin(X))
                miny = min(miny, np.amin(Y))

                maxx = max(maxx, np.amax(X))
                maxy = max(maxy, np.amax(Y))

                points = ax[i,j].scatter(X, Y, c=colors, s=40, alpha=0.6)
                allPltAnnot.append((points, htmlData))

        combinedPoints = []
        combinedHTML = []

        for tooltipdata in allPltAnnot:
            combinedPoints.append(tooltipdata[0])
            combinedHTML.append(tooltipdata[1])

        # remove tick labels
        #for axi in ax.flat:
        #    for axis in [axi.xaxis, axi.yaxis]:
        #        axis.set_major_formatter(plt.NullFormatter())

        minx = math.floor(minx)
        miny = math.floor(miny)

        maxx = math.ceil(maxx)
        maxy = math.ceil(maxy)

        stepSizeX = int((maxx-minx)/10.0)
        stepSizeY = int((maxy-miny)/10.0)

        for axi in ax.flat:

            axi.set_ylim([miny, maxy])
            axi.set_xlim([minx, maxx])

            axi.xaxis.set_ticks([x for x in range(minx, maxx, stepSizeX)])
            axi.yaxis.set_ticks([x for x in range(miny, maxy, stepSizeY)])

            alllabels = axi.xaxis.get_ticklabels()

            for tick in axi.xaxis.get_major_ticks():
                tick.label._text = str(tick._loc)
                tick.label.set_rotation(65)

        # Here we connect the linked brush plugin
        tooltip = MultiAxesPointHTMLTooltip( combinedPoints, combinedHTML, voffset=10, hoffset=10, css=tooltipCSS)
        mpld3.plugins.connect(fig, tooltip)
        mpld3.plugins.connect(fig, mpld3.plugins.MousePosition(fontsize=14))

        mpld3.plugins.connect(fig, mpld3.plugins.LinkedBrush(points))
        mpld3.show()

        return

    def get_id(self, obj, suffix="", prefix="el", warn_on_invalid=True):
        """Get a unique id for the object"""
        if not suffix:
            suffix = ""
        if not prefix:
            prefix = ""
        import os

        objid = prefix + str(os.getpid()) + str(id(obj)) + suffix

        if warn_on_invalid and not self.html_id_ok(objid):
            print('"{0}" is not a valid html ID. This may cause problems')

        return objid

    def html_id_ok(self, objid, html5=False):
        """Check whether objid is valid as an HTML id attribute.

        If html5 == True, then use the more liberal html5 rules.
        """
        import re

        if html5:
            return not re.search('\s', objid)
        else:
            return bool(re.match("^[a-zA-Z][a-zA-Z0-9\-\.\:\_]*$", objid))









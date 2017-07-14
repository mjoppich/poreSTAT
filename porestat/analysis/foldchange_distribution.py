import argparse
from collections import defaultdict

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

from scipy.optimize import curve_fit
from scipy.misc import factorial
class FoldChangeDistributionFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(FoldChangeDistributionFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('foldchange', help='expls help')
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

        return FoldChangeAnalysis(simArgs)



class FoldChangeAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(FoldChangeAnalysis, self).__init__(args)
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


        dataShape = cntData.shape


        allFCs = defaultdict(list)

        for i in range(dataShape[0]):
            for j in range(dataShape[0]):

                if (i==j):
                    continue


                X = cntData[j, :]
                Y = cntData[i, :]


                for k in range(0, len(unionGenes)):

                    if i != j and X[k] == Y[k] and X[k] > 0:
                        print(unionGenes[k])

                    if (X[k] == 0) or (Y[k] == 0):
                        continue

                    #print(Y[k]/X[k])
                    #print(math.log(Y[k] / X[k], 2))
                    allFCs[ (dataRows[i], dataRows[j]) ].append( math.log(Y[k]/X[k], 2) )

                print(allFCs[ (dataRows[i], dataRows[j]) ])

                entries, bin_edges, patches = plt.hist(allFCs[ (dataRows[i], dataRows[j]) ], bins=50, normed=True)

                # calculate binmiddles
                bin_middles = 0.5 * (bin_edges[1:] + bin_edges[:-1])

                # poisson function, parameter lamb is the fit parameter
                def poisson(k, lamb):
                    return (lamb ** k / factorial(k)) * np.exp(-lamb)

                # fit with curve_fit
                parameters, cov_matrix = curve_fit(poisson, bin_middles, entries)

                # plot poisson-deviation with fitted parameter
                x_plot = np.linspace(-100, 20, 100)

                plt.plot(x_plot, poisson(x_plot, *parameters), 'r-', lw=2)
                plt.show()

                plt.show()

        return










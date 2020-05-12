import argparse
from collections import defaultdict

import HTSeq
import matplotlib

from porestat.plots.poreplot import PorePlot, MultiAxesPointHTMLTooltip

from porestat.plots.plotconfig import PlotConfig, PlotSaveTYPE
from porestat.utils.EnrichmentDF import EnrichmentDF
from ..utils.DataFrame import DataFrame, DataRow, ExportTYPE
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory, PSToolException, PSToolInterface

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Files import fileExists, eprint
from scipy import stats
import sys, math
import numpy as np
from matplotlib import pyplot as plt
import random
import os

from scipy.optimize import curve_fit


class FoldChangeSimilarityFactory(PSToolInterfaceFactory):
    def __init__(self, parser, subparsers, which):

        super(FoldChangeSimilarityFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-d', '--diffreg', nargs='+', type=str, help='poreSTAT diffreg results',
                            required=True)
        parser.add_argument('-m', '--methods', type=str, nargs='+', default=['edgeR', 'DESeq'])

        parser = PlotConfig.addParserArgs(parser)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return FoldChangeSimilarity(simArgs)


class FoldChangeSimilarity(ParallelPSTInterface):
    def __init__(self, args):

        super(FoldChangeSimilarity, self).__init__(args)

        self.counts = None
        self.condData = EnrichmentDF()

    def _makePropDict(self):

        return None

    def readCounts(self, args):

        for countsFile in args.counts:
            if not fileExists(countsFile):
                raise PSToolException("Read info file does not exist: " + str(countsFile))

        counts = {}
        for countsFile in args.counts:
            counts[countsFile] = DataFrame.parseFromFile(countsFile)

        return counts

    def readDiffRegs(self, args):

        for countsFile in args.diffreg:
            if not fileExists(countsFile):
                raise PSToolException("Diffreg file does not exist: " + str(countsFile))

        for diffFile in args.diffreg:
            df = EnrichmentDF.parseFromFile(diffFile)

    def prepareInputs(self, args):
        return []

    def execParallel(self, data, environment):

        return None

    def joinParallel(self, existResult, newResult, oEnvironment):

        return None

    def makeResults(self, parallelResult, oEnvironment, args):

        allDiffRegData = {}
        allDiffRegSims = defaultdict(dict)
        conditions = set()

        for file in args.diffreg:

            thisData = EnrichmentDF(DataFrame.parseFromFile(file))
            condPair = tuple(thisData.getConditions())

            for cond in condPair:
                conditions.add(cond)

            allDiffRegData[condPair] = thisData

            for method in args.methods:

                methodFCs = []
                for x in thisData.getColumn(method + "_log2FC"):
                    if x != None and x!= 'None':
                        methodFCs.append(abs(float(x)))

                average = sum(methodFCs) / len(methodFCs)

                allDiffRegSims[method][condPair] = average


        allConditions = sorted(list(conditions))

        for method in allDiffRegSims:

            sims = np.zeros( (len(allConditions), len(allConditions)) )

            for condPair in allDiffRegSims[method]:
                sims[ allConditions.index(condPair[0]), allConditions.index(condPair[1]) ] = allDiffRegSims[method][condPair]
                sims[allConditions.index(condPair[1]), allConditions.index(condPair[0])] = allDiffRegSims[method][
                    condPair]

            PorePlot.heat_map_cluster(sims, allConditions, allConditions, "Similarity: " + str(method), "", pltcfg=args.pltcfg)





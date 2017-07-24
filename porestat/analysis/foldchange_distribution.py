import argparse
from collections import defaultdict

import HTSeq
import matplotlib

from porestat.plots.poreplot import PorePlot, MultiAxesPointHTMLTooltip

from porestat.plots.plotconfig import PlotConfig, PlotSaveTYPE
from porestat.utils.EnrichmentDF import EnrichmentDF
from ..utils.DataFrame import DataFrame, DataRow, ExportTYPE
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Files import fileExists
from scipy import stats
import sys, math
import numpy as np
from matplotlib import pyplot as plt
import random
import os


from scipy.optimize import curve_fit
from scipy.misc import factorial
class FoldChangeDistributionFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(FoldChangeDistributionFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('foldchange', help='expls help')
        parser.add_argument('-c', '--counts', nargs='+', type=str, default=None, help='counts summary file', required=False)
        parser.add_argument('-d', '--diffreg', nargs='+', type=str, default=None, help='poreSTAT diffreg results', required=False)
        parser.add_argument('-v', '--no-analysis', dest='noanalysis', action='store_true', default=False)
        parser.add_argument('-m', '--methods', type=str, nargs='+', default=['edgeR', 'DESeq'])

        parser.add_argument('-o', '--output', type=str, help='output location, default: std out', default=sys.stdout)
        parser.add_argument('-r', '--rscript', type=str, help='path to Rscript', default='/usr/bin/Rscript')

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


        if not args.counts == None and not args.noanalysis:

            counts = self.readCounts(args)

            vConds = sorted([x for x in counts])

            createdComparisons = defaultdict(list)
            conditions = []

            for valueSource in ['coverage', 'read_counts']:
                self.condData = EnrichmentDF()

                for condition in vConds:

                    condData = counts[condition]

                    geneNames = condData.getColumnIndex('gene')
                    geneCounts = condData.getColumnIndex(valueSource)

                    condRow = condData.toDataRow(geneNames,geneCounts)

                    condition = condition.split("/")[-2]
                    conditions.append(condition)

                    self.condData.addCondition(condRow, condition)

                print("Running for conditions: " + str(vConds))

                createdComparisons[valueSource] += self.condData.runDEanalysis( args.output, prefix = valueSource, rscriptPath=args.rscript, methods=args.methods )

            self.prepareHTMLOut(createdComparisons, args)


        if args.diffreg != None:

            createdComparisons = defaultdict(list)
            conditions = set()

            for file in args.diffreg:

                df = EnrichmentDF.parseFromFile(file)
                valueSource = self.getValueSource(df)

                conditions.add(df.getHeader()[2])
                conditions.add(df.getHeader()[3])

                createdComparisons[valueSource].append(file)

            self.prepareHTMLOut(createdComparisons, conditions, args)

    def getValueSource(self, df):
        return df.data[0][0]

    def prepareHTMLOut(self, createdComparisons, conditions, args):
        
        for valueSource in createdComparisons:
            self.condData.printResult(args.output, outputPrefix=valueSource, conditions=conditions, files=createdComparisons[valueSource])










import argparse
from collections import defaultdict

import HTSeq
import matplotlib
from porestat.utils.ArgParseExt import FolderType

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
from scipy.misc import factorial
class FoldChangeDistributionFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(FoldChangeDistributionFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), default=None, help='counts summary file', required=False)
        parser.add_argument('-d', '--diffreg', nargs='+', type=str, default=None, help='poreSTAT diffreg results', required=False)
        parser.add_argument('-v', '--no-analysis', dest='noanalysis', action='store_true', default=False)
        parser.add_argument('-m', '--methods', type=str, nargs='+', default=['NOISeq', 'DESeq'])

        parser.add_argument('-o', '--output', type=FolderType('w'), help='output location, default: std out', required=True)
        parser.add_argument('-r', '--rscript', type=argparse.FileType('r'), help='path to Rscript', default='/usr/bin/Rscript')

        parser = PlotConfig.addParserArgs(parser)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        if (not PSToolInterface.hasArgument('counts', args) and not PSToolInterface.hasArgument('diffreg', args) ) or (args.counts == None and args.diffreg == None):
            raise argparse.ArgumentParser().error("error: Either counts [--counts] or diffreg results [--diffreg] must be set!")

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
            if not fileExists(countsFile.name):
                raise PSToolException("Read info file does not exist: " + str(countsFile))

        counts = {}
        for countsFile in args.counts:
            counts[countsFile.name] = DataFrame.parseFromFile(countsFile.name)

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

                createdComparisons[valueSource] += self.condData.runDEanalysis( args.output, prefix = valueSource, rscriptPath=args.rscript.name, methods=args.methods )

            self.prepareHTMLOut(createdComparisons, conditions, args)


        if args.diffreg != None:

            createdComparisons = defaultdict(list)
            conditions = set()

            for file in args.diffreg:

                df = EnrichmentDF.parseFromFile(file)
                valueSource = self.getValueSource(df)

                conditions += df.getConditions()

                createdComparisons[valueSource].append(file)

            self.prepareHTMLOut(createdComparisons, conditions, args)

    def getValueSource(self, df):
        return df.data[0][1]

    def prepareHTMLOut(self, createdComparisons, conditions, args):
        
        for valueSource in createdComparisons:
            self.condData.printResult(args.output, prefix=valueSource, conditions=conditions, files=createdComparisons[valueSource])












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


class FoldChangeFeatureCountsDistributionFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):
        super(FoldChangeFeatureCountsDistributionFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which + ' help')
        parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), default=None,
                            help='counts summary file', required=False)
        parser.add_argument('-p', '--prefixes', nargs='+', type=str, default=None,
                            help='counts summary file', required=False)

        parser.add_argument('-co', '--conditions', nargs='+', type=str, action='append', default=None,
                            help='counts summary file', required=False)

        parser.add_argument('-v', '--no-analysis', dest='noanalysis', action='store_true', default=False)
        parser.add_argument('-m', '--methods', type=str, nargs='+', default=['DESeq', 'edgeR', 'limma', 'msEmpiRe'])

        parser.add_argument('-o', '--output', type=FolderType('w'), help='output location, default: std out',
                            required=True)
        parser.add_argument('-r', '--rscript', type=argparse.FileType('r'), help='path to Rscript',
                            default='/usr/bin/Rscript')

        parser.add_argument('-e', '--enhanced', type=argparse.FileType('r'), default=None)

        parser = PlotConfig.addParserArgs(parser)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        if (not PSToolInterface.hasArgument('counts', args) and not PSToolInterface.hasArgument('diffreg', args)) or (
                args.counts == None and args.diffreg == None):
            raise argparse.ArgumentParser().error(
                "error: Either counts [--counts] or diffreg results [--diffreg] must be set!")

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return FoldChangeFeatureCountsAnalysis(simArgs)


class FoldChangeFeatureCountsAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(FoldChangeFeatureCountsAnalysis, self).__init__(args)

        self.counts = None
        self.condData = EnrichmentDF()

    def _makePropDict(self):

        return None

    def readCounts(self, args):

        featureCountsColumns = ["Geneid",	"Chr",	"Start",	"End",	"Strand",	"Length"]

        counts = defaultdict(lambda: list())

        condition2samples = defaultdict(list)

        for idx, countFile in enumerate(args.counts):

            countFilePrefix = args.prefixes[idx]
            df = DataFrame.parseFromFile(countFile.name, skipChar='#')

            allheaders = df.getHeader()
            sampleHeaders = [x for x in allheaders if not x in featureCountsColumns]

            for sample in sampleHeaders:
                condition2samples[sample].append(countFilePrefix + sample)


            for condGroup in args.conditions:
                condName = condGroup[0]
                for condElement in condGroup:

                    print(condName, condElement)

                    subDf = df.selectColumns({"Geneid": "gene", condElement: "count"})
                    subDf.setFilepath(os.path.abspath(condElement))

                    counts[condName].append(subDf)

        return counts, condition2samples

    def readDiffRegs(self, args):

        # TODO what did I want to do with this argument?
        for diffFile in args.diffreg:
            df = EnrichmentDF.parseFromFile(diffFile)

    def prepareInputs(self, args):
        return []

    def execParallel(self, data, environment):

        return None

    def joinParallel(self, existResult, newResult, oEnvironment):

        return None

    def loadEnhancement(self, fileE):

        if fileE == None:
            print("Not loading gene name enhancements")
            return {}

        print("Loading gene name enhancements", fileE.name)

        ens2sym = {}
        for line in fileE:
            line = line.strip().split("\t")

            if not line[0].startswith("ENS"):
                continue

            ensemblID = line[0]
            geneSymbol = line[1]

            if len(geneSymbol) == 0:
                continue

            ens2sym[ensemblID] = geneSymbol

        return ens2sym


    def makeResults(self, parallelResult, oEnvironment, args):

        if not args.counts == None:

            """
            counts is a defaultdict(list) for each condition name with maybe multiple samples
            """

            geneEnhancement = self.loadEnhancement(args.enhanced)

            counts, cond2samples = self.readCounts(args)

            vConds = sorted([x for x in counts])

            createdComparisons = defaultdict(list)
            conditions = []

            for valueSource in ['count']:
                self.condData = EnrichmentDF()
                replicates = {}

                for condition in vConds:

                    condData = counts[condition]

                    condReplicates = []
                    for condDataSample in condData:
                        geneNames = condDataSample.getColumnIndex('gene')
                        geneCounts = condDataSample.getColumnIndex(valueSource)

                        condRow = condDataSample.toDataRow(geneNames, geneCounts)

                        sampleName = condDataSample.filepath
                        conditions.append(sampleName)

                        condReplicates.append(sampleName)

                        self.condData.addCondition(condRow, sampleName)

                    replicates[condition] = condReplicates

                print("Running for conditions: " + str(vConds))

                createdComparisons[valueSource] += self.condData.runDEanalysis(args.output, prefix=valueSource,
                                                                               rscriptPath=args.rscript.name,
                                                                               methods=args.methods,
                                                                               replicates=replicates,
                                                                               noDErun=args.noanalysis,
                                                                               enhanceSymbol=geneEnhancement)

            self.prepareHTMLOut(createdComparisons, replicates, args)


    def getValueSource(self, df):
        return df.data[0][1]

    def prepareHTMLOut(self, createdComparisons, replicates, args):

        for valueSource in createdComparisons:

            allComparisons = createdComparisons[valueSource]

            condPair2File = {}
            for x in allComparisons:
                condPair2File[(x[0], x[1])] = x[2]

            print("Comparisons")
            for x in condPair2File:
                print(x, condPair2File[x])

            self.condData.printResult(args.output, prefix=valueSource, conditionPair2File=condPair2File,
                                      replicates=replicates)  # conditions=conditions, files=createdComparisons[valueSource])

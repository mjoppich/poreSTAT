import argparse
import HTSeq
import matplotlib

from porestat.plots.poreplot import PorePlot, MultiAxesPointHTMLTooltip

from porestat.plots.plotconfig import PlotConfig
from ..utils.DataFrame import DataFrame, DataRow
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter, defaultdict
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
        parser.add_argument('-d', '--counts', nargs='+', type=str, help='counts summary file', required=False)
        parser.add_argument('-o', '--output', type=str, help='output location, default: std out', default=None)

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

    def prepareInputs(self, args):
        return []

    def execParallel(self, data, environment):

        return None


    def joinParallel(self, existResult, newResult, oEnvironment):

        return None


    def getMeasuredEvidences(self, df):


        evidences = set()
        measures = df.getHeader()[1:]

        for i in range(0, len(measures), 2):
            evidence = measures[i]
            evidences.add(evidence)

        return evidences

    def makeResults(self, parallelResult, oEnvironment, args):

        counts = self.readCounts(args)

        vConds = sorted([x for x in counts])

        createdComparisons = defaultdict(list)
        conditions = []
        evidences = set()

        for condition in vConds:
            df = counts[condition]

            dfMeasures = self.getMeasuredEvidences(df)
            evidences = evidences.union(dfMeasures)

            condition = condition.split("/")[-2]
            conditions.append(condition)

        for evidence in evidences:

            similarities = {}
            simData = np.zeros((len(conditions), len(conditions)))

            for i in range(0, len(vConds)):

                cond1File = vConds[i]
                cond1 = conditions[i]
                cond1Data = self.getConditionData(counts[cond1File], evidence)

                for j in range(i+1, len(conditions)):
                    cond2File = vConds[j]

                    cond2 = conditions[j]
                    cond2Data = self.getConditionData(counts[cond2File], evidence)

                    similarity = self.calcSimilarity(cond1Data, cond2Data)
                    similarities[ (cond1, cond2) ] = similarity
                    simData[ i,j ] = similarity
                    simData[j,i] = similarity


            for condPair in similarities:
                print(str(condPair) + " " + str(similarities[condPair]))

            PorePlot.heat_map_cluster(simData, conditions, conditions, "Similarity", "cond sim")

    def getConditionData(self, df, evidence):

        geneNames = df.getColumnIndex('gene')
        geneCounts = df.getColumnIndex(evidence)
        condRow = df.toDataRow(geneNames, geneCounts)

        return condRow

    def calcSimilarity(self, cond1, cond2):

        unionGenes = set()

        cond1Header = cond1.getHeader()
        cond2Header = cond2.getHeader()

        for x in cond1Header:

            if x in cond2Header:
                unionGenes.add(x)

        totalCond1 = sum([cond1[x] for x in cond1Header])
        totalCond2 = sum([cond2[x] for x in cond2Header])

        commonPairs = []
        dist = 0

        for ugene in unionGenes:

            cond1val = cond1[ugene] / totalCond1
            cond2val = cond2[ugene] / totalCond2

            commonPairs = (cond1val, cond2val)

            dist = math.sqrt(sum([math.pow((commonPairs[0]-commonPairs[1]), 2)]))
            #dist += commonPairs[0] * commonPairs[1]

        return 1.0-dist







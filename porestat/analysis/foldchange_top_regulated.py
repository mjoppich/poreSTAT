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
from scipy.misc import factorial


class FoldChangeTopRegulatedFactory(PSToolInterfaceFactory):
    def __init__(self, parser, subparsers, which):

        super(FoldChangeTopRegulatedFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-d', '--diffreg', nargs='+', type=str, help='poreSTAT diffreg results',
                            required=True)
        parser.add_argument('-m', '--methods', type=str, nargs='+', default=['edgeR', 'DESeq'])
        parser.add_argument('-n', '--top', type=int, default=20)
        parser.add_argument('-o', '--output', type=str, required=True)

        parser = PlotConfig.addParserArgs(parser)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return FoldChangeTopRegulatedAnalysis(simArgs)


class FoldChangeTopRegulatedAnalysis(ParallelPSTInterface):
    def __init__(self, args):

        super(FoldChangeTopRegulatedAnalysis, self).__init__(args)

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

        def parseNones(row):

            ret = [None] * len(row)
            for i in range(0, len(row)):
                if row[i] != 'None':
                    ret[i] = row[i]

            return ret

        topGenes = Counter()

        for file in args.diffreg:

            thisData = EnrichmentDF(DataFrame.parseFromFile(file))
            thisData.applyToRow( parseNones )

            condPair = tuple(thisData.getConditions())

            for cond in condPair:
                conditions.add(cond)

            allDiffRegData[condPair] = thisData

            for method in args.methods:

                methodFCs = []

                pvals = thisData.toDataRow(thisData.getColumnIndex('id'), thisData.getColumnIndex(method + "_RAW.PVA"))
                genepval = [(x[0], float(x[1])) for x in pvals.to_pairs() if x[1] != None]

                genepval.sort(key=lambda x: x[1])

                for i in range(0, args.top):
                    topGenes[ genepval[i][0] ] += 1

        outputDF = DataFrame()

        geneIDidx = outputDF.addColumn('gene_id')
        countIdx = outputDF.addColumn('count')
        linkIdx = outputDF.addColumn('link')


        for (gene, count) in topGenes.most_common():

            geneRow = DataRow.fromDict( {
                'gene_id': gene,
                'count': count,
                'link': "<a href='http://www.uniprot.org/uniprot/?query="+gene+"&sort=score' target='_blank'>UniProt</a>",
            })

            outputDF.addRow(geneRow)

        outputDF.export(args.output, ExportTYPE.HTML)








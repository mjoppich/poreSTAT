import argparse
from collections import defaultdict

import HTSeq
import matplotlib

from porestat.plots.poreplot import PorePlot, MultiAxesPointHTMLTooltip

from porestat.plots.plotconfig import PlotConfig
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
import mpld3
import random
import os


class EnrichmentDF(DataFrame):

    def __init__(self):
        super(EnrichmentDF, self).__init__()

        self.idCol = self.addColumn('id')

    def addCondition(self, condData, condName):

        isDataRow = isinstance(condData, DataRow)
        isDict = isinstance(condData, dict)

        if not (isDataRow or isDict):
            raise ValueError("parameter not of type datarow or dict")

        if isDataRow:
            return self._addConditionFromDataRow(condData, condName)

        if isDict:
            return self._addConditionFromDict(condData, condName)


    def _addConditionFromDict(self, condData, condName):
        dr = DataRow.fromDict(condData)
        return self._addConditionFromDataRow(dr, condName)


    def _addConditionFromDataRow(self, condData, condName):

        newColIdx = self.addColumn(condName, None)

        setFoundKeys = set()

        for i in range(0, len(self.data)):
            row = self.data[i]
            rowID = row[ self.idCol ]

            if rowID in condData.getHeader():
                lrow = list(row)
                lrow[newColIdx] = condData[rowID]
                self.data[i] = tuple(lrow)
                setFoundKeys.add(rowID)

        for geneID in condData.getHeader():
            if not geneID in setFoundKeys:

                lrow = [None] * len(self.column2idx)
                lrow[0] = geneID
                lrow[newColIdx] = condData[geneID]

                self.data.append(tuple(lrow))


    def writeEnrichmentBrowserFiles(self, prepData, exprFile, pdataFile, fdataFile):

        self.writeExpressionFile(exprFile, prepData)
        self.writepDataFile(pdataFile, prepData)
        self.writefDataFile(fdataFile, prepData)



    def prepareEBData(self, cond1, cond2):

        header = ['gene', cond1, cond2]
        prepData = [tuple(header)]

        cond1Idx = self.getColumnIndex(cond1)
        cond2Idx = self.getColumnIndex(cond2)

        for row in self:

            count1 = int(row[cond1Idx])
            count2 = int(row[cond2Idx])

            gene = row[0]

            if count1 == None or count2 == None:
                continue

            prepData.append( (gene, count1, count2) )

        return prepData


    def writepDataFile(self, file, prepData):

        with open(file, 'w') as pdataFile:

            conds = prepData[0][1:]
            cnt = 0
            for cond in conds:
                pdataFile.write(cond + "\t" + str(cnt) + "\n")
                cnt += 1


    def writefDataFile(self, file, prepData):

        with open(file, 'w') as fdataFile:

            for line in prepData[1:]:
                fdataFile.write( line[0] + "\t" + line[0] + "\n")

    def writeExpressionFile(self, file, prepData):

        with open(file, 'w') as exprFile:

            for line in prepData[1:]:
                exprFile.write( "\t".join([str(x) for x in line[1:]]) + "\n")



    def runDEanalysis(self, conditions=None):

        if conditions == None:
            conditions = self.getHeader()[1:]

        base = "/tmp/eb/"

        exprFile = base + "expr"
        pdataFile = base + "p_data"
        fdataFile = base + "f_data"
        outFileBase = base + "out_data"

        scriptPath = os.path.dirname(os.path.abspath(__file__)) + "/../data/de_rseq.R"


        for i in range(0, len(conditions)):
            cond1 = conditions[i]
            for j in range(i+1, len(conditions)):
                cond2 = conditions[j]

                prepData = self.prepareEBData(cond1, cond2)
                self.writeEnrichmentBrowserFiles(prepData, exprFile, pdataFile, fdataFile)

                condResult = {}

                for method in ['limma', 'edgeR', 'DESeq']:
                    outFile = outFileBase + "_" + method

                    os.system("/home/proj/biosoft/software/R/R-3.2.2/bin/Rscript "+scriptPath+" "+exprFile+" "+pdataFile+" "+fdataFile+" "+method+" " + outFile)

                    methDF = DataFrame.parseFromFile(outFile)
                    condResult[method] = methDF

                geneNames = self.getColumnIndex('id')
                cond1Counts = self.toDataRow(geneNames, self.getColumnIndex(cond1))
                cond2Counts = self.toDataRow(geneNames, self.getColumnIndex(cond2))


                compDF = EnrichmentDF()
                compDF.addCondition(cond1Counts, cond1)
                compDF.addCondition(cond2Counts, cond2)

                usedMethod = []
                for method in condResult:

                    methDF = condResult[method]

                    if methDF == None:
                        continue

                    l2FCTitle = method + "_log2FC"
                    rawpTitle = method + "_RAW.PVA"
                    adjpTitle = method + "_ADJ.PVAL"

                    l2FCdata = methDF.toDataRow( methDF.getColumnIndex('GENE.ID'), methDF.getColumnIndex('log2FC') )
                    rawPdata = methDF.toDataRow(methDF.getColumnIndex('GENE.ID'), methDF.getColumnIndex('RAW.PVAL'))
                    adjPdata = methDF.toDataRow(methDF.getColumnIndex('GENE.ID'), methDF.getColumnIndex('ADJ.PVAL'))

                    compDF.addCondition(l2FCdata, l2FCTitle)
                    compDF.addCondition(rawPdata, rawpTitle)
                    compDF.addCondition(adjPdata, adjpTitle)

                print(compDF)

                addInfo = compDF.addColumn("EnsemblBacteria")

                def addInfoFunc(x):
                    gene = x[0]
                    link = "<a target='_blank' href='http://bacteria.ensembl.org/Helicobacter_pylori_p12/Gene/Summary?g="+gene+"'>EnsemblBacteria</a>"

                    if not gene.startswith("HP"):
                        link = ""

                    x[addInfo] = link
                    return tuple(x)

                compDF.applyToRow( addInfoFunc )

                (head, body) = compDF.export("/tmp/eb/table.html", ExportTYPE.HTML_STRING)


                return


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
        self.condData = EnrichmentDF()

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

        for condition in self.counts:

            condData = self.counts[condition]

            geneNames = condData.getColumnIndex('gene')
            geneCounts = condData.getColumnIndex('coverage')

            condRow = condData.toDataRow(geneNames,geneCounts)

            condition = condition.split("/")[8]

            self.condData.addCondition(condRow, condition)

        self.condData.runDEanalysis()










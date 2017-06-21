from collections import OrderedDict

import math

import time

from porestat.plots.poreplot import PorePlot

from porestat.plots.plotconfig import PlotConfig

from .ParallelPTTInterface import ParallelPSTReportableInterface
from .PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

import matplotlib.pyplot as plt

import argparse
class QualityPositionFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(QualityPositionFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('qual_pos', help='expls help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-p', '--no-plot', action='store_true', default=False)
        parser.add_argument('-u', '--user_run', dest='groupByUser', action='store_true', default=False)
        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return QualityPosition(simArgs)

class QualityPosition(ParallelPSTReportableInterface):

    def __init__(self, args):

        super(QualityPosition, self).__init__( args )
        self.qualTypes = [chr(x) for x in range(33, 126+1)]

    def _makePropDict(self):

        propDict = {}
        propDict['QUALS'] = {}
        propDict['QUAL_SIMPLE'] = Counter()
        propDict['USER_RUN_NAME'] = set()
        propDict['READ_COUNT'] = 0

        return propDict

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def handleEntity(self, fileObj, localEnv, globalEnv):

        runid = fileObj.runID()

        if not runid in localEnv:
            localEnv[runid] = self._makePropDict()

        propDict = localEnv[runid]
        propDict['READ_COUNT'] += 1
        propDict['USER_RUN_NAME'].add(fileObj.user_filename_input())

        fastq = fileObj.getFastQ()

        if fastq != None:

            qualDict = propDict['QUALS']

            for i in range(0, len(fastq.qual)):

                found_qual = fastq.qual[i]

                if not i in qualDict:
                    qualDict[i] = Counter()

                qualDict[i][found_qual] += 1

                propDict['QUAL_SIMPLE'][found_qual] += 1

            propDict['QUALS'] = qualDict

        return localEnv


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'FILES']

        qualObservations = []

        for x in self.qualTypes:
            qualObservations.append(x)

        qualObservations.append("NTs")

        for x in self.qualTypes:
            qualObservations.append(x + "%")

        allobservations = {}
        for runid in parallelResult:

            props = parallelResult[runid]

            run_user_name = props['USER_RUN_NAME']
            fileCount = props['READ_COUNT']
            simpleQualCounter = props['QUAL_SIMPLE']
            qualCounter = props['QUALS']

            observations = {

                'RUNID': runid,
                'USER_RUN_NAME': ",".join(run_user_name),
                'FILES': fileCount,
                'QUALPOS': qualCounter,
                'QUAL_SIMPLE': simpleQualCounter
            }

            key = ",".join(run_user_name) if self.hasArgument('groupByRunName', args) and args.groupByRunName else runid

            if key in allobservations:
                allobservations[key] = mergeDicts(allobservations[key], observations)
            else:
                allobservations[key] = observations

        sortedruns = sorted([x for x in allobservations])

        print("\t".join(makeObservations + qualObservations))

        plotData = OrderedDict()

        for runid in sortedruns:

            allobs = []
            for x in makeObservations:
                allobs.append( str(allobservations[runid][x]) )

            simpleQualObs = allobservations[runid]['QUAL_SIMPLE']

            totalCount = 0
            for qual in self.qualTypes:
                count = simpleQualObs[qual]
                totalCount += count
                allobs.append( str(count) )

            allobs.append(str(totalCount))

            for qual in self.qualTypes:
                count = simpleQualObs[qual] / totalCount
                allobs.append( "{0:.4g}".format(count) )

            print("\t".join(allobs))

            # make plot for runid

            if not self.hasArgument('no_plot', args) or args.no_plot == False:
                # pos -> qual -> count
                qualCounter = allobservations[runid]['QUALPOS']
                plotData[runid] = qualCounter

        if not self.hasArgument('no_plot', args) or args.no_plot == False:
            self.plotQualCounter(plotData, args)


    def plotQualCounter(self, qualCounters, args):


        foundLengths = set()

        for runid in qualCounters:
            for length in qualCounters[runid]:
                foundLengths.add(length)


        steps = 10

        foundLengths = sorted(foundLengths)
        maxLength = max(foundLengths)
        step = maxLength / steps;

        allPlotData = {}
        stepLabels = []

        minQual = None
        maxQual = None

        for runid in qualCounters:

            allDataPlot = []
            qualCounter = qualCounters[runid]

            for i in range(0, steps):

                stepMin = i*step
                stepMax = stepMin + step

                stepLabels.append( "{0:.0f}-{1:.0f}".format(stepMin, stepMax) )

            # create empty counter for each step
            for i in range(0, steps):
                allDataPlot.append( Counter() )


            for pos in qualCounter:

                bin = divmod(pos, math.ceil(step))

                for q in qualCounter[pos]:
                    allDataPlot[bin[0]][q] += qualCounter[pos][q]

            explodedData = []
            for counterData in allDataPlot:

                dataVal = []
                for q in counterData:
                    dataVal = dataVal + [ord(q)] * counterData[q]

                    minQual = ord(q) if (minQual == None) or (minQual > ord(q)) else minQual
                    maxQual = ord(q) if (maxQual == None) or (maxQual < ord(q)) else maxQual

                explodedData.append(dataVal)
            
            allPlotData[runid] = explodedData

        minQual = 32
        maxQual = 127

        def axManipulation(axis):
            axis.axes.get_xaxis().set_ticks( [i for i in range(1, len(stepLabels)+1)] )
            axis.axes.get_xaxis().set_ticklabels( stepLabels, rotation=90 )

            axis.axes.get_yaxis().set_ticks([i for i in range(minQual, maxQual+1)])
            axis.axes.get_yaxis().set_ticklabels([str(chr(i)) for i in range(minQual, maxQual+1)])


        allAxisManip = [axManipulation] * len(allPlotData)

        PorePlot.plotViolin( allPlotData, None, "Quality Position Distribution", axisManipulation=allAxisManip, pltcfg=args.pltcfg )






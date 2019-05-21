import sys
from collections import OrderedDict, defaultdict

import math

import time

from porestat.plots.poreplot import PorePlot, PlotDirectionTYPE

from porestat.plots.plotconfig import PlotConfig, PlotSaveTYPE

from .ParallelPTTInterface import ParallelPSTReportableInterface
from .PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

import matplotlib.pyplot as plt

import argparse
class QualityPositionFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(QualityPositionFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-p', '--no-plot', action='store_true', default=False)
        parser.add_argument('-u', '--user_run', dest='user_run', action='store_true', default=False)
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

            key = self.makeKey(run_user_name, args, runid)

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
        maxLength = max(foundLengths)+1
        step = maxLength / steps;

        allPlotData = {}
        stepLabels = []

        minQual = None
        maxQual = None

        for i in range(0, steps):
            stepMin = i * step
            stepMax = stepMin + step

            stepLabels.append("{0:.0f}-{1:.0f}".format(stepMin, stepMax))

        allRunIDs = []

        for runid in qualCounters:

            allDataPlot = []
            qualCounter = qualCounters[runid]

            # create empty counter for each step
            for i in range(0, steps):
                allDataPlot.append( Counter() )


            for pos in qualCounter:

                bin = divmod(pos, math.ceil(step))

                if (bin[0] > len(allDataPlot)):
                    bin = (len(allDataPlot)-1, 0)

                for q in qualCounter[pos]:
                    allDataPlot[bin[0]][q] += qualCounter[pos][q]

            explodedData = []
            for counterData in allDataPlot:

                elemCounts = [counterData[q] for q in counterData if counterData[q] > 20]

                if len(elemCounts) > 0:
                    scaleFactor = min(elemCounts)
                else:
                    scaleFactor = 1

                print("min values", scaleFactor, [(ord(q), q) for q in counterData], [int(counterData[q]/scaleFactor) for q in counterData])


                dataVal = []
                for q in counterData:

                    cdQ = counterData[q]
                    dataVal = dataVal + [ord(q)] * int(cdQ/scaleFactor)

                    minQual = ord(q) if (minQual == None) or (minQual > ord(q)) else minQual
                    maxQual = ord(q) if (maxQual == None) or (maxQual < ord(q)) else maxQual

                explodedData.append(sorted(dataVal))

                if len(dataVal) > 0:
                    print(min(dataVal), max(dataVal))
            
            allPlotData[runid] = explodedData
            allRunIDs.append(runid)

        minQual = 32
        maxQual = 127

        def axManipulation(axis):
            axis.axes.get_xaxis().set_ticks( [i for i in range(1, len(stepLabels)+1)] )
            axis.axes.get_xaxis().set_ticklabels( stepLabels, rotation=90 )

            axis.axes.get_yaxis().set_ticks([i for i in range(minQual, maxQual+1)])
            axis.axes.get_yaxis().set_ticklabels([str(chr(i)) for i in range(minQual, maxQual+1)])


        allAxisManip = [axManipulation] * len(allPlotData)

        #PorePlot.plotViolin( allPlotData, None, "Quality Position Distribution", axisManipulation=allAxisManip, pltcfg=args.pltcfg )
        self.makePlot( allPlotData, stepLabels, sorted(allRunIDs), args.pltcfg )

    def makePlot(self, allPlotData, steps, runIDs, pltcfg):

        lengthsByStep = defaultdict(lambda: defaultdict(list))

        for runid in runIDs:

            for i in range(0, len(allPlotData[runid])):
                lengthsByStep[ runid ][ steps[i] ] = allPlotData[runid][i]


        for runid in runIDs:

            for x in lengthsByStep[runid]:
                maxval = min(100, len(lengthsByStep[runid][x]))
                print(x, len(lengthsByStep[runid][x]), lengthsByStep[runid][x][0:maxval])

            PorePlot.plotViolinSNS(lengthsByStep[runid], None, 'Quality at Position ('+runid+")", pltcfg=pltcfg,
                                plotDirection=PlotDirectionTYPE.VERTICAL, xTitle="Read Length Range", yTitle="Read Quality Scores")

        """

        pltcfg.startPlot()

        figSize = (10, (2+len(runIDs)) * len(steps))

        fig, ax = plt.subplots(nrows=len(steps), ncols=1, sharex=True, sharey=True, figsize=figSize)

        for i in range(0, len(steps)):

            step = steps[i]
            plotData = lengthsByStep[step]

            dataToPlot = []
            for runid in runIDs:
                if runid in plotData:
                    dataToPlot.append(plotData[runid])
                else:
                    dataToPlot.append([])

            print("Plotting step: " + str(step))
            PorePlot.plotSingleViolin(dataToPlot, "Range " + str(step) + "bp", ax[i], vert=False)
            ax[i].set_yticks([x for x in range(1, len(runIDs)+1)])
            ax[i].set_yticklabels(runIDs)
            ax[i].yaxis.tick_right()

            for x in ax[i].get_yticklabels():
                print(x)


        pltcfg.makePlot(figHeight=None)

        """






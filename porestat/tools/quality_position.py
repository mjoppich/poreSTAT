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
        parser.add_argument('-p', '--no-plot', nargs='?', action='store_true', default=False)
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

            key = ",".join(run_user_name) if args.groupByUser else runid

            if key in allobservations:
                allobservations[key] = mergeDicts(allobservations[key], observations)
            else:
                allobservations[key] = observations

        sortedruns = sorted([x for x in allobservations])

        print("\t".join(makeObservations + qualObservations))

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

            if self.hasArgument('no_plot', args) and args.no_plot == False:
                # pos -> qual -> count
                qualCounter = observations['QUALPOS']
                self.plotQualCounter(qualCounter, args)


    def plotQualCounter(self, qualCounter, args):


        foundLengths = set()
        for length in qualCounter:
            foundLengths.add(length)

        foundLengths = sorted(foundLengths)
        maxLength = max(foundLengths)
        step = maxLength / 10;

        allDataPlot = []
        stepLabels = []

        minQual = None
        maxQual = None

        for i in range(0, 10):

            stepMin = i*step
            stepMax = stepMin + step

            stepLabels.append( "{0:.0f}-{1:.0f}".format(stepMin, stepMax) )

            allData = Counter()
            for pos in qualCounter:

                if stepMin <= pos and pos < stepMax:

                    for k in qualCounter[pos]:
                        allData[k] += qualCounter[pos][k]

            dataVal = []
            dataPos = []

            for x in allData:
                dataVal = dataVal + [ord(x)] * allData[x]

                minQual = ord(x) if (minQual == None) or (minQual > ord(x)) else minQual
                maxQual = ord(x) if (maxQual == None) or (maxQual < ord(x)) else maxQual


            allDataPlot.append(dataVal)

        minQual = 32
        maxQual = 127

        fig, ax = plt.subplots()
        ax.violinplot(allDataPlot, showmeans=True, showextrema=True, showmedians=True)

        ax.axes.get_xaxis().set_ticks( [i for i in range(1, len(stepLabels)+1)] )
        ax.axes.get_xaxis().set_ticklabels( stepLabels, rotation=90 )

        ax.axes.get_yaxis().set_ticks([i for i in range(minQual, maxQual+1)])
        ax.axes.get_yaxis().set_ticklabels([str(chr(i)) for i in range(minQual, maxQual+1)])

        args.pltcfg.makePlot()






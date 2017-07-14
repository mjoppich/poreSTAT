from collections import defaultdict

from ..plots.plotconfig import PlotConfig
from ..plots.poreplot import PorePlot, PlotDirectionTYPE

from .ParallelPTTInterface import ParallelPSTReportableInterface
from .PTToolInterface import PSToolInterfaceFactory, PSToolException
from ..utils.Stats import calcN50

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

class LengthHistogramFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(LengthHistogramFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('hist', help='expls help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-p', '--plot', nargs='?', type=bool, const=True, default=False, help='issue plot?', required=False)
        parser.add_argument('-u', '--user-run', dest='user_run', action='store_true', default=False)
        parser.add_argument('-q', '--read-type', dest='read_type', action='store_true', default=False, help='add type subplots')


        parser.add_argument('-c', '--combined', dest='combineRuns', action='store_true', default=False)
        parser.add_argument('-v', '--violin', dest='violin', action='store_true', default=False)

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return LengthHistogram(simArgs)

class LengthHistogram(ParallelPSTReportableInterface):

    def __init__(self, args):

        super(LengthHistogram, self).__init__( args )

    def _makePropDict(self):

        propDict = {}
        propDict['USER_RUN_NAME'] = set()
        propDict['READ_COUNT'] = 0
        propDict['LENGTHS'] = []

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
            propDict['LENGTHS'].append((len(fastq), fileObj.type))

        return localEnv

    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        if parallelResult == None:
            raise PSToolException('No valid result generated.')

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'FILES', 'TOTAL_LENGTH', 'N50', 'L50']

        allobservations = {}
        for runid in parallelResult:

            props = parallelResult[runid]

            run_user_name = props['USER_RUN_NAME']
            fileCount = props['READ_COUNT']
            lengthObversations = props['LENGTHS']

            onlyLengths = [x[0] for x in lengthObversations]

            (n50, l50) = calcN50(onlyLengths)

            if not args.read_type:
                lengthObversations = onlyLengths

            observations = {

                'RUNID': runid,
                'USER_RUN_NAME': ",".join(run_user_name),
                'FILES': fileCount,
                'TOTAL_LENGTH': sum(onlyLengths),
                'N50': n50,
                'L50': l50,
                'LENGTHS': lengthObversations
            }

            key = ",".join(run_user_name) if self.hasArgument('user_run', args) and args.user_run else runid

            if key in allobservations:
                allobservations[key] = mergeDicts(allobservations[key], observations)
            else:
                allobservations[key] = observations

        sortedruns = sorted([x for x in allobservations])

        print("\t".join(makeObservations))

        plotData = {}

        if self.hasArgument('combineRuns', args) and args.combineRuns:

            lengthData = []
            labels = []

            for runid in sortedruns:
                lengthData += allobservations[runid]['LENGTHS']
                labels.append(runid)

                self.printObservation(makeObservations, allobservations[runid])

            plotLabel = ",".join(labels)
            plotData[ plotLabel ] = lengthData

        else:

            for runid in sortedruns:

                self.printObservation(makeObservations, allobservations[runid])
                plotData[runid] =  allobservations[runid]['LENGTHS']

        if self.hasArgument('read_type', args) and args.read_type:

            newPlotData = {}

            for runid in plotData:
                lengthsByType = defaultdict(list)

                for x in plotData[runid]:
                    lengthsByType[x[1]].append(x[0])

                for readtype in lengthsByType:
                    newid = runid + "_" + readtype
                    newPlotData[newid] = lengthsByType[readtype]

            plotData = newPlotData

        if self.hasArgument('violin', args) and not args.violin:
            PorePlot.plotHistogram(plotData, None, 'Length Histogram for ', xlabel="Read Length", ylabel="Read Count", pltcfg=args.pltcfg)
        else:
            PorePlot.plotViolin(plotData, None, 'Length Histogram', pltcfg=args.pltcfg, plotDirection=PlotDirectionTYPE.HORIZONTAL)



    def printObservation(self, makeObservations, thisObservation):

        allobs = []
        for x in makeObservations:
            allobs.append(str(thisObservation[x]))

        print("\t".join(allobs))


from ..plots.poreplot import PorePlot
from ..plots.plotconfig import PlotConfig

from .ParallelPTTInterface import ParallelPSTReportableInterface
from .PTToolInterface import PSToolInterfaceFactory
from ..utils.Stats import calcN50

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

class YieldPlotFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(YieldPlotFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-mr', '--mreads', nargs='+', type=str, help='multi-read files', required=False)

        parser.add_argument('-u', '--user-run', action='store_true', default=False, help='groups all runs by user run field')

        parser.add_argument('--no-read-type-subplot', action='store_true', default=False, help='do not add type subplots')
        parser.add_argument('--separate-subplots', action='store_true', default=False,
                            help='add subplots by read type, but make separate figure')

        parser.add_argument('-rc', '--read-count', dest='read_count', action='store_true', default=False)

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return YieldPlot(simArgs)


class YieldData:

    def __init__(self, runid):
        self.runid = runid
        self.byType = {}
        self.alldata = None

    def addData(self, data, dataType):

        self.byType[dataType] = data

    def getData(self, dataType = None):

        if dataType == None:

            if self.alldata != None:
                return self.alldata
            else:
                return []

        else:

            return self.byType[dataType]




class YieldPlot(ParallelPSTReportableInterface):

    def __init__(self, args):

        super(YieldPlot, self).__init__( args )

    def _makePropDict(self):

        propDict = {}
        propDict['USER_RUN_NAME'] = set()
        propDict['READ_COUNT'] = 0
        propDict['TIME_LENGTHS'] = []

        return propDict

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)


    def handleEntity(self, entity, localEnv, globalEnv):
        runid = entity.runID()

        if not runid in localEnv:
            localEnv[runid] = self._makePropDict()

        propDict = localEnv[runid]
        propDict['READ_COUNT'] += 1
        propDict['USER_RUN_NAME'].add(entity.user_filename_input())

        fastq = entity.getFastQ()

        if fastq != None:
            timeOfCreation = entity.readCreateTime() - entity.getExperimentStartTime()
            readLength = len(fastq)
            readType = entity.type

            propDict['TIME_LENGTHS'].append((timeOfCreation, readLength, readType))

        return localEnv

    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        if parallelResult == None:
            return

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'FILES', 'TOTAL_LENGTH', 'N50', 'L50']

        allobservations = {}

        for runid in parallelResult:

            props = parallelResult[runid]

            run_user_name = props['USER_RUN_NAME']
            fileCount = props['READ_COUNT']

            if self.hasArgument('read_count', args) and args.read_count:
                lengthObversations = [1 for x in props['TIME_LENGTHS']] # sums number of reads up
            else:
                lengthObversations = [x[1] for x in props['TIME_LENGTHS']]

            (n50, l50) = calcN50(lengthObversations)

            observations = {

                'RUNID': runid,
                'USER_RUN_NAME': ",".join(run_user_name),
                'FILES': fileCount,
                'TOTAL_LENGTH': sum(lengthObversations),
                'N50': n50,
                'L50': l50,
                'TIME_LENGTHS': props['TIME_LENGTHS']
            }

            key = self.makeKey(run_user_name, args, runid)

            if key in allobservations:
                allobservations[key] = mergeDicts(allobservations[key], observations)
            else:
                allobservations[key] = observations

        sortedruns = sorted([x for x in allobservations])

        print("\t".join(makeObservations))

        for runid in sortedruns:

            allobs = []
            for x in makeObservations:
                allobs.append(str(allobservations[runid][x]))

            print("\t".join(allobs))

        self.makePlot( allobservations, args )


    def makePlot(self, data, args):


        # prepare data
        makeReadTypeData = (not (self.hasArgument('no_read_type_subplot', args)) or not args.no_read_type_subplot)
        rawPlotData = self.prepareData(data, makeReadTypeData)

        if self.hasArgument('read_count', args) and args.read_count:
            yDescription = "Cumulative Read Count"
        else:
            yDescription = "Cumulative BP"




        plotData = {}

        for runid in rawPlotData:
            plotData[runid] = rawPlotData[runid].getData()

        if args.separate_subplots:
            PorePlot.yieldPlot(plotData, "Yield Plot (all reads)", "Time since exp. start", yDescription, pltcfg=args.pltcfg)
            plotData = {}

        for type in Fast5TYPE:
            for runid in rawPlotData:
                plotData[runid + " " + str(type)] = rawPlotData[runid].getData(dataType=type)

            if args.separate_subplots:

                totalElements = 0
                for x in plotData:
                    totalElements += len(plotData[x])

                if totalElements == 0:
                    plotData = {}
                    continue

                PorePlot.yieldPlot(plotData, "Yield Plot (" + str(type) + ")", "Time since exp. start", yDescription, pltcfg=args.pltcfg)
                plotData = {}

        if not args.separate_subplots:
            PorePlot.yieldPlot(plotData, "Yield Plot", "Time since exp. start", yDescription, pltcfg=args.pltcfg)

    def prepareData(self, data, makeReadTypeData):

        timeLengthData = {}

        for runid in data:

            locTL = sorted(data[runid]['TIME_LENGTHS'], key=lambda x: x[0])

            cumL = 0
            cumLocTL = []
            for tpl in locTL:
                locT = tpl[0]
                locL = tpl[1]

                cumLocTL.append((locT, locL + cumL))

                cumL += locL

            runData = YieldData(runid)
            runData.alldata = cumLocTL

            """

            BY READ TYPE

            """

            if makeReadTypeData:

                cumLocTLbyType = {}
                cumLByType = {}

                locTL = sorted(data[runid]['TIME_LENGTHS'], key=lambda x: x[0])

                for tpl in locTL:

                    locT = tpl[0]
                    locL = tpl[1]
                    locType = tpl[2]

                    if not locType in cumLByType:
                        cumLByType[locType] = 0
                        cumLocTLbyType[locType] = [(0, 0)]

                    cumL = cumLByType[locType]
                    cumLocTL = cumLocTLbyType[locType]

                    cumL += locL
                    cumLocTL.append((locT, cumL))

                    cumLByType[locType] = cumL
                    cumLocTLbyType[locType] = cumLocTL

                for type in Fast5TYPE:
                    cumData = cumLocTLbyType[type] if type in cumLocTLbyType else []
                    runData.addData(cumData, type)

            timeLengthData[runid] = runData

        return timeLengthData
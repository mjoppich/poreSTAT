from ..plots.poreplot import PorePlot, TimestampTimeFormatter

from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory
from ..utils.Stats import calcN50

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

class YieldPlotFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(YieldPlotFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('yield', help='expls help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-p', '--plot', nargs='?', type=bool, const=True, default=False, help='issue plot?', required=False)
        parser.add_argument('-u', '--user_run', dest='groupByRunName', action='store_true', default=False)
        parser.add_argument('-q', '--read_type', dest='addTypeSubplot', action='store_true', default=False, help='add type subplots')
        parser.add_argument('-v', '--violin', dest='violin', action='store_true', default=False)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return YieldPlot(simArgs)

class YieldPlot(ParallelPSTInterface):

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

    def execParallel(self, data, environment):

        counterRunID = {}

        f5folder = Fast5Directory(data)

        iFilesInFolder = 0

        for file in f5folder.collect():

            runid = file.runID()

            iFilesInFolder += 1

            if not runid in counterRunID:
                counterRunID[runid] = self._makePropDict()

            propDict = counterRunID[runid]
            propDict['READ_COUNT'] += 1
            propDict['USER_RUN_NAME'].add( file.user_filename_input() )

            fastq = file.getFastQ()

            if fastq != None:

                timeOfCreation = file.readCreateTime() - file.getExperimentStartTime()
                readLength = len(fastq)
                readType = file.type

                propDict['TIME_LENGTHS'].append( (timeOfCreation, readLength, readType) )

        print("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return counterRunID


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'FILES', 'TOTAL_LENGTH', 'N50', 'L50']

        allobservations = {}
        for runid in parallelResult:

            props = parallelResult[runid]

            run_user_name = props['USER_RUN_NAME']
            fileCount = props['READ_COUNT']
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

            key = ",".join(run_user_name) if args.groupByRunName else runid

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

        timeAndLength = []
        labels = []

        timeLengthData = {}

        for runid in data:

            locTL = sorted(data[runid]['TIME_LENGTHS'], key=lambda tup: tup[0])

            cumL = 0

            cumLocTL = []
            for tpl in locTL:

                locT = tpl[0]
                locL = tpl[1]
                cumLocTL.append( (locT, locL+cumL) )

                cumL += locL

            locTL = cumLocTL

            timeLengthData[runid] = locTL

            if args.addTypeSubplot:

                cumLocTLbyType = {}
                cumLByType = {}

                for tpl in locTL:

                    locT = tpl[0]
                    locL = tpl[1]
                    locType = tpl[2]

                    if not locType in cumLByType:
                        cumLByType[locType] = 0
                        cumLocTLbyType[locType] = []


                    cumL = cumLByType[locType]
                    cumLocTL = cumLocTLbyType[locType]

                    cumLocTL.append( (locT, locL+cumL) )
                    cumL += locL

                    cumLByType[locType] = cumL
                    cumLocTLbyType[locType] = cumLocTL

                for dataType in cumLocTLbyType:

                    label = runid + "_" + str(dataType)

                    timeLengthData[ label ] = cumLocTLbyType[dataType]

        PorePlot.yieldPlot(timeLengthData, "Yield Plot", "Cumulative BP", "Time")



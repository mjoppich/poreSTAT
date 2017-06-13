import argparse
import os

import sys

from numpy import genfromtxt
from porestat.plots.plotconfig import PlotConfig

from ..utils.Utils import mergeDicts
from ..utils.Stats import calcN50

from .ParallelPTTInterface import ParallelPSTReportableInterface
from .PTToolInterface import PSToolInterfaceFactory, PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from ..plots.poreplot import PorePlot

class ChannelOccupancyFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(ChannelOccupancyFactory, self).__init__(parser, self._addParser(subparsers))

    def _addParser(self, subparsers):
        parser = subparsers.add_parser('occ', help='occ help')

        parser = subparsers.add_parser('occ', help='occ help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan')
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-e', '--experiments', nargs='+', type=str, help='experiments to list')
        parser.add_argument('-u', '--user_run', dest='groupByRunName', action='store_true', default=False)
        parser.add_argument('-t', '--tsv', nargs='?', action='store', type=argparse.FileType('w'), const=sys.stdout, default=None)

        parser.add_argument('-q', '--read_type', nargs='+', type=str, choices=[x for x in Fast5TYPE.str2type], help='read types ('+ ",".join([x for x in Fast5TYPE.str2type]) +')')

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return ChannelOccupancy(simArgs)



class ChannelOccupancy(ParallelPSTReportableInterface):

    def __init__(self, args):
        super(ChannelOccupancy, self).__init__(args)

        this_dir, this_filename = os.path.split(__file__)
        chipLayout = genfromtxt(this_dir + '/../data/chip_layout.csv', delimiter=',', dtype=int)

        self.channelCount = chipLayout.shape[0] * chipLayout.shape[1]


    def _makePropDict(self):

        propDict = {}
        propDict['USER_RUN_NAME'] = set()
        propDict['READ_COUNT'] = 0
        propDict['CHANNELS'] = {}

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

        channelID = fileObj.channelID()
        channelDict = propDict['CHANNELS']
        timeOfCreation = fileObj.readCreateTime() - fileObj.getExperimentStartTime()

        fastq = fileObj.getFastQ()

        readLength = 0
        if fastq != None:
            readLength = len(fastq)

        readType = str(fileObj.type)

        if globalEnv.read_type != None and not readType in globalEnv.read_type:
            return localEnv

        if channelID in channelDict:
            channelDict[channelID].append((timeOfCreation, readLength))
        else:
            channelDict[channelID] = [(timeOfCreation, readLength)]

    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult

    def prepareEnvironment(self, args):

        return args


    def makeResults(self, parallelResult, oEnvironment, args):

        if parallelResult == None:
            raise PSToolException("Something went wrong. Empty result.")

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'FILES', 'TOTAL_LENGTH', 'N50', 'L50']

        allobservations = {}
        for runid in parallelResult:

            props = parallelResult[runid]

            run_user_name = props['USER_RUN_NAME']
            fileCount = props['READ_COUNT']

            lengthObversations = []
            for channelID in props['CHANNELS']:
                lengthObversations = lengthObversations + [x[1] for x in props['CHANNELS'][channelID]]

            (n50, l50) = calcN50(lengthObversations)

            observations = {

                'RUNID': runid,
                'USER_RUN_NAME': ",".join(run_user_name),
                'FILES': fileCount,
                'TOTAL_LENGTH': sum(lengthObversations),
                'N50': n50,
                'L50': l50,
                'CHANNELS': props['CHANNELS']
            }

            key = ",".join(run_user_name) if self.hasArgument('groupByRunName', args) and args.groupByRunName else runid

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

        printTsv = self.hasArgument('tsv', args) and args.tsv

        if printTsv:
            # print header

            vChannels = [ x for x in range(1, self.channelCount+1) ]
            args.tsv.write( 'run_name' + "\t" + "\t".join([str(x) for x in vChannels]) + "\n" )

        for runid in sortedruns:

            if printTsv:
                self.printChannelHistogram(args.tsv, runid, vChannels, allobservations[runid]['CHANNELS'])
            else:
                self.makePlot(runid, allobservations[runid]['CHANNELS'], args)

        if printTsv:
            args.tsv.flush()
            args.tsv.close()

    def printChannelHistogram(self, file, runid, vChannels, channelDict):

        channel2rl = {}

        for channelID in channelDict:
            channel2rl[ channelID ] = [str(x[1]) for x in channelDict[channelID]]

        file.write(runid + "\t")

        channelVec = []
        for channelID in vChannels:

            if not channelID in channel2rl:
                channelVec.append( "" )
            else:
                channelVec.append( ",".join(channel2rl[channelID]) )

        file.write("\t".join(channelVec) + "\n")


    def makePlot(self, runKey, channelDict, args):

        channel2rl = {}

        for channelID in channelDict:
            channel2rl[channelID] = [x[1] for x in channelDict[channelID]]

        PorePlot.plotLoadOut(channel2rl, runKey, pltcfg=args.pltcfg)
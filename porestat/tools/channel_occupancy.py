import argparse
import os

import sys
from collections import Counter
from collections import defaultdict

from numpy import genfromtxt
from ..plots.plotconfig import PlotConfig

from ..utils.Utils import mergeDicts
from ..utils.Stats import calcN50

from .ParallelPTTInterface import ParallelPSTReportableInterface
from .PTToolInterface import PSToolInterfaceFactory, PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE, Fast5TYPEAction
from ..plots.poreplot import PorePlot
from ..utils.DataFrame import DataFrame, ExportTYPEAction, ExportTYPE, DataRow

class ChannelOccupancyFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(ChannelOccupancyFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan')
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-e', '--experiments', nargs='+', type=str, help='experiments to list')
        parser.add_argument('-u', '--user_run', dest='user_run', action='store_true', default=False)
        parser.add_argument('--print-histogram', action='store_true', default=False)
        
        parser.add_argument('-o', '--output', nargs='?', type=str, default=None, const=None)
        parser.add_argument('-ot', '--output-type', nargs='?', action=ExportTYPEAction, default=ExportTYPE.CSV)

        parser.add_argument('-q', '--read-type', dest='read_type', action=Fast5TYPEAction)

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
        self.allChannels = [i for i in range(1, self.channelCount+1)]

    def _makePropDict(self):

        propDict = {}
        propDict['USER_RUN_NAME'] = set()
        propDict['READ_COUNT'] = 0
        propDict['CHANNELS'] = defaultdict(list)

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

            key = self.makeKey(run_user_name, args, runid)

            if key in allobservations:
                allobservations[key] = mergeDicts(allobservations[key], observations)
            else:
                allobservations[key] = observations

        sortedruns = sorted([x for x in allobservations])

        self.printStats(makeObservations, allobservations)

        for runid in sortedruns:
            self.makePlot(runid, allobservations[runid]['CHANNELS'], args)

        if self.hasArgument('print-histogram', args) and args.print_histogram == True:
            self.printChannelHistogram(args, allobservations)



    def printStats(self, makeObservations, allobservations):
        print("\t".join(makeObservations))

        for runid in sorted([x for x in allobservations]):

            allobs = []
            for x in makeObservations:
                allobs.append(str(allobservations[runid][x]))

            print("\t".join(allobs))

    def printChannelHistogram(self, args, allObservations):

        vHeaders = ['run_name'] + self.allChannels

        outputFrame = DataFrame()
        outputFrame.addColumns(vHeaders)

        for runID in allObservations:

            channelDict = allObservations[runID]

            channel2rl = defaultdict(list)
            channel2rl['run_name'] = runID

            allChannels = channelDict['CHANNELS']

            for channelID in self.allChannels:
                channel2rl[ channelID ] = [x[1] for x in allChannels[channelID]]

            runData = DataRow.fromDict(channel2rl)
            outputFrame.addRow(runData)

        outputFrame.export(args.output, args.output_type)


    def makePlot(self, runKey, channelDict, args):

        channel2rl = defaultdict(list)

        for channelID in channelDict:
            channel2rl[channelID] = [x[1] for x in channelDict[channelID]]

        PorePlot.plotLoadOut(channel2rl, runKey, pltcfg=args.pltcfg)
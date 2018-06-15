from ..plots.plotconfig import PlotConfig
from ..plots.poreplot import PorePlot

from ..utils.Utils import mergeCounter
from ..utils.Files import eprint

from .ParallelPTTInterface import ParallelPSTInterface
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE

import argparse

from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

from ..utils.PickleArgparse import FileType

class Environment(object):
    pass

class TimelineReadsFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(TimelineReadsFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-o', '--out', action='store', type=argparse.FileType('w'), default=None)

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        simArgs = self._makeArguments(args)

        return TimelineReads(simArgs)

class TimelineReads(ParallelPSTInterface):

    def __init__(self, args):

        super(TimelineReads, self).__init__( args )


    def prepareEnvironment(self, args):

        oEnv = Environment()
        return oEnv

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def execParallel(self, data, environment):

        f5folder = Fast5Directory(data)

        iFilesInFolder = 0

        createdReadsTime = Counter()
        createdReadsBasecalledTime = Counter()

        for file in f5folder.collect():

            iFilesInFolder += 1

            winnerseq = file.getFastQ()
            creationTime = file.readCreateTime()

            if winnerseq != None:
                createdReadsBasecalledTime[creationTime] += 1
            else:
                createdReadsTime[creationTime] += 1

        eprint("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return [createdReadsTime, createdReadsBasecalledTime]

    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = [Counter(), Counter()]


        existResult = [
            mergeCounter(existResult[0], newResult[0]),
            mergeCounter(existResult[1], newResult[1])
        ]

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        readsTime = parallelResult[0]
        readsBCTime = parallelResult[1]

        PorePlot.plotTimeLine([readsBCTime, readsTime], ['basecalled', 'non-basecalled'], 'Timeline plot')








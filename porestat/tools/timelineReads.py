
from ..plots.poreplot import PorePlot
from ..utils.Utils import eprint
from .ParallelPTTInterface import ParallelPTTInterface
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE

import argparse

from collections import Counter


class Environment(object):
    pass



class TimelineReads(ParallelPTTInterface):

    def __init__(self, parser, subparsers):

        super(TimelineReads, self).__init__(parser, self.__addParser(subparsers))


    def __addParser(self, subparsers):

        parser_expls = subparsers.add_parser('time', help='timeline help')
        parser_expls.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser_expls.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser_expls.add_argument('-o', '--out', action='store', type=argparse.FileType('w'), default=None)

        parser_expls.set_defaults(func=self.exec)

        return parser_expls

    def prepareEnvironment(self, args):

        oEnv = Environment()
        return oEnv

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def execParallel(self, procID, environment, data):

        f5folder = Fast5Directory(data)

        iFilesInFolder = 0

        createdReadsTime = Counter()
        createdReadsBasecalledTime = Counter()

        for file in f5folder.collect():

            runid = file.runID()

            iFilesInFolder += 1

            runid = file.runID()

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

        PorePlot.plotTimeLine([readsBCTime, readsTime], ['basecalled', 'non-basecalled'])








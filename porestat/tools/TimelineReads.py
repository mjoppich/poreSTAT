from ..utils.Utils import eprint
from .ParallelPTTInterface import ParallelPTTInterface
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE

from collections import Counter

import numpy as np
from matplotlib.mlab import csv2rec
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
from matplotlib.ticker import Formatter
import sys, argparse, os

class Environment(object):
    pass

class MyFormatter(Formatter):
    def __init__(self, dates, fmt='%Y-%m-%d %H:%M:%S'):
        self.dates = dates
        self.fmt = fmt

    def __call__(self, x, pos=0):
        'Return the label for time x at position pos'
        ind = int(np.round(x))
        if ind >= len(self.dates) or ind < 0:
            return ''

        return self.dates[ind].strftime(self.fmt)

class TimelineReads(ParallelPTTInterface):

    def __init__(self, parser, subparsers):

        super(TimelineReads, self).__init__(parser, self.__addParser(subparsers))


    def __addParser(self, subparsers):

        parser_expls = subparsers.add_parser('seq', help='seq help')
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

            winnerseq = file.winner()
            creationTime = file.createTime()

            if winnerseq != None:
                createdReadsBasecalledTime[creationTime] += 1
            else:
                createdReadsTime[creationTime] += 1

        eprint("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return [createdReadsTime, createdReadsBasecalledTime]

    def mergeCounter(self, counter1, counter2):

        mergedCounter = Counter()

        for x in counter1:
            mergedCounter[x] = counter1[x]

        for x in counter2:
            mergedCounter[x] += counter2[x]

        return mergedCounter


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = [Counter(), Counter()]


        existResult = [
            self.mergeCounter(existResult[0], newResult[0]),
            self.mergeCounter(existResult[1], newResult[1])
        ]

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        readsTime = parallelResult[0]
        readsBCTime = parallelResult[1]

        timePoints = readsTime.keys().sort()
        timePointsBC = readsBCTime.keys().sort()

        valPoints = [readsTime[x] for x in timePoints]
        valPointsBC = [readsBCTime[x] for x in timePointsBC]

        fig, ax = plt.subplots()

        formatter = MyFormatter(timePoints)
        ax.xaxis.set_major_formatter(formatter)
        ax.plot(np.arange(len(timePoints)), valPoints, 'o-')

        formatter = MyFormatter(timePoints)
        ax.xaxis.set_major_formatter(formatter)
        ax.plot(np.arange(len(timePointsBC)), valPointsBC, 'x-')

        fig.autofmt_xdate()
        plt.show()








import argparse

import sys

from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE, Fast5TYPEAction


from collections import OrderedDict

class ReadInfoFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(ReadInfoFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('info', help='expls help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)

        parser.add_argument('-q', '--read-type', nargs='+', dest='read_type', action=Fast5TYPEAction, default=None)
        parser.add_argument('-u', '--user-run', dest='user_run', action='store_true', default=False)
        parser.add_argument('-e', '--experiments', nargs='+', type=str,
                            help='run ids of experiments to be extracted. if --user_run, give user_run_name s',
                            required=False)

        def fileOpener( filename ):
            open(filename, 'w').close()
            return filename

        parser.add_argument('-o', '--output', type=fileOpener, help='output location, default: std out', default=sys.stdout)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return ReadInfo(simArgs)

import os

class Environment(object):
    pass


class ReadInfo(ParallelPSTInterface):

    def __init__(self, args):

        super(ReadInfo, self).__init__(args)

        dReadSummary = OrderedDict()

        dReadSummary['READ_ID'] = lambda file: file.readID()
        dReadSummary['READ_NAME'] = lambda file: file.sequenceName()

        dReadSummary['CHANNEL_ID'] = lambda file: file.channelID()
        dReadSummary['READ_NUMBER'] = lambda file: file.readNumber()

        dReadSummary['TYPE'] = lambda file: file.type
        dReadSummary['READ_LENGTH'] = lambda file: file.sequenceLength()
        dReadSummary['AVG_QUALITY'] = lambda file: "0"
        dReadSummary['TIME'] = lambda file: file.readCreateTime()

        dReadSummary['USER_RUN_NAME'] = lambda file: file.user_filename_input()
        dReadSummary['RUN_ID'] = lambda file: file.runID()

        self.endl = os.linesep
        self.dReadSummary = dReadSummary
        self.dObservations = [x for x in self.dReadSummary]


    def _makePropDict(self):
        return None

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def execParallel(self, data, environment):

        foundReads = []

        f5folder = Fast5Directory(data)

        iFilesInFolder = 0

        for file in f5folder.collect():

            iFilesInFolder += 1

            dReadInfo = self.makeReadInfo(file)

            foundReads.append(dReadInfo)

        print("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return foundReads

    def makeReadInfo(self, readFile):

        dReadInfo = {}

        for key in self.dReadSummary:

            dReadInfo[key] = str( self.dReadSummary[key](readFile) )

        return dReadInfo

    def prepareEnvironment(self, args):

        env = Environment()
        env.output = args.output

        header = "\t".join(self.dObservations) + self.endl

        self.writeLinesToOutput(env.output, [header], mode='w')

        return env


    def joinParallel(self, existResult, newResult, environment):

        if newResult == None:
            return None

        allLines = []

        for readEntry in newResult:

            readinfo = [ readEntry[x] for x in self.dObservations ]
            toWrite = "\t".join(readinfo) + self.endl
            allLines.append(toWrite)

        if not allLines is None and len(allLines) > 0:
            self.writeLinesToOutput(environment.output, allLines)

        return None



    def makeResults(self, parallelResult, env, args):
        self.closeOutput(env.output)

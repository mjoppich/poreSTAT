from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

class StatsSummaryFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(StatsSummaryFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser_expls = subparsers.add_parser('summary', help='read performance summary')
        parser_expls.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser_expls.add_argument('-r', '--reads', type=str, help='minion read folder', required=False)
        parser_expls.set_defaults(func=self._prepObj)

        return parser_expls

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return StatsSummary(simArgs)

class StatsSummary(ParallelPSTInterface):

    def __init__(self, args):

        super(StatsSummary, self).__init__( args )

    def exec(self, args):

        folders = self.manage_folders_reads(args)

        counterRunID = {}

        for folder in folders:

            f5folder = Fast5Directory(folder)

            for file in f5folder.collect():

                runid = file.runID()

                if not runid in counterRunID:
                    counterRunID[runid] = 1
                else:
                    counterRunID[runid] += 1

        print(counterRunID)

        print("executed LS")
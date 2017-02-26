from .PTToolInterface import PTToolInterface
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory

class StatsSummary(PTToolInterface):

    def __init__(self, parser, subparsers):

        super(StatsSummary, self).__init__(parser, self.__addParser(subparsers))

    def __addParser(self, subparsers):

        parser_expls = subparsers.add_parser('stats', help='expls help')
        parser_expls.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser_expls.add_argument('-r', '--reads', type=str, help='minion read folder', required=False)
        parser_expls.set_defaults(func=self.exec)

        return parser_expls


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
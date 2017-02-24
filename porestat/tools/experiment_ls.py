from .PTToolInterface import PTToolInterface
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory

class Experiment_ls(PTToolInterface):

    def __init__(self, parser, subparsers):

        super(Experiment_ls, self).__init__(parser, self.__addParser(subparsers))

    def __addParser(self, subparsers):

        parser_expls = subparsers.add_parser('expls', help='expls help')
        parser_expls.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser_expls.add_argument('-r', '--reads', type=str, help='minion read folder', required=False)
        parser_expls.set_defaults(func=self.exec)

        return parser_expls


    def exec(self, args):

        if (args.folders == None and args.reads == None):
            print("error: Either folders or reads must be set!")
            self.subparser.print_help()
            exit(-1)


        if args.folders != None:
            folders = args.folders
        else:
            folders = []

        if args.reads != None:
            folders.append( '' )

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
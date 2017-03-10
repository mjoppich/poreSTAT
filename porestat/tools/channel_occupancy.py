from .PTToolInterface import PTToolInterface
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from porestat.plots.poreplot import PorePlot

class Channel_occupancy(PTToolInterface):

    def __init__(self, parser, subparsers):

        super(Channel_occupancy, self).__init__(parser, self.__addParser(subparsers))


    def __addParser(self, subparsers):

        parser_chocc = subparsers.add_parser('occ', help='occ help')
        parser_chocc.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan')
        parser_chocc.add_argument('-r', '--reads', type=str, help='minion read folder', required=False)
        parser_chocc.add_argument('-e', '--experiments', nargs='+', type=str, help='experiments to list')
        parser_chocc.set_defaults(func=self.exec)

        return parser_chocc


    def exec(self, args):


        folders = self.manage_folders_reads(args)

        experiments = None
        if not args.experiments == None:
            experiments = set(args.experiments)

        exp2cl = {}
        runid2filename = {}

        for folder in folders:

            print("Checking folder: " + folder)

            f5folder = Fast5Directory(folder)

            for file in f5folder.collect():

                if file.type == Fast5TYPE.UNKNOWN:
                    continue

                runid = file.runID()
                runid2filename[runid] = file.user_filename_input()

                if (experiments != None and len(experiments) > 0 and runid in experiments) or (experiments == None):

                    if not runid in exp2cl:
                        exp2cl[runid] = {}

                    channelLengths = exp2cl[runid]

                    channelID = file.channelID()
                    readLength = len(file.winnerFQ())

                    if not channelID in channelLengths:
                        channelLengths[channelID] = [readLength]
                    else:
                        channelLengths[channelID].append(readLength)

                    exp2cl[runid] = channelLengths

        for runid in exp2cl:

            myset = set()

            for x in exp2cl[runid]:
                myset.add(x)

            print(runid)
            print(runid2filename[runid])
            print(myset)
            print(len(myset))

            PorePlot.plotLoadOut(exp2cl[runid], pores=(32,16))
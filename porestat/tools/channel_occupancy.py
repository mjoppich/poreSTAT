from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from porestat.plots.poreplot import PorePlot

class ChannelOccupancyFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(ChannelOccupancyFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('occ', help='occ help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan')
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-e', '--experiments', nargs='+', type=str, help='experiments to list')
        parser.add_argument('-t', '--by_type', dest='byType', action='store_true', default=False)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return ChannelOccupancy(simArgs)



class ChannelOccupancy(ParallelPSTInterface):

    def __init__(self, args):
        super(ChannelOccupancy, self).__init__(args)

    def exec(self):


        folders = self.manage_folders_reads(self.args)

        experiments = None
        if not self.args.experiments == None:
            experiments = set(self.args.experiments)

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
                    readLength = len(file.getFastQ())

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

            PorePlot.plotLoadOut(exp2cl[runid], runid2filename[runid])
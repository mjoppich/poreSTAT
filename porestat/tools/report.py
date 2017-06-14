import os

from porestat.tools.length_histogram import LengthHistogram
from porestat.tools.nucleotide_distribution import NucleotideDistribution
from porestat.tools.quality_distribution import QualityDistribution
from porestat.utils.Files import makePath
from .channel_occupancy import ChannelOccupancy
from .yield_plot import YieldPlot


from ..plots.plotconfig import PlotConfig
from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

from collections import OrderedDict


class ReportFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(ReportFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('summary', help='read performance summary')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)

        parser.add_argument('-o', '--output', type=str, help='output folder for report', required=True)
        parser.add_argument('-n', '--output-name', type=str, help='output name', required=False)

        parser.add_argument('--no-read-type-subplot', dest='addTypeSubplot', action='store_false', default=True, help='do not add type subplots')
        parser.add_argument('-q', '--read-type', nargs='+', type=str, choices=[x.value for x in Fast5TYPE], help='read types ('+ ",".join([x.value for x in Fast5TYPE]) +')')
        parser.add_argument('-u', '--user-run', dest='groupByRunName', action='store_true', default=False)


        parser = PlotConfig.addParserArgs(parser)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return ReportAnalysis(simArgs)


class ReportAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(ReportAnalysis, self).__init__( args )

        self.dReporters = OrderedDict([
            ('OCCUPANCY', ChannelOccupancy(args)),
            ('YIELD', YieldPlot(args)),
            ('LENGTH', LengthHistogram(args)),
            ('NUCLEOTIDES', NucleotideDistribution(args)),
            ('QUALITY DISTRIBUTION', QualityDistribution(args)),
            ('QUALITY BY POSITION', QualityDistribution(args)),
        ])

        print("Output folder: " + str(args.output))
        print("Output name:   " + str(args.output_name))

        if args.output_name == None:
            args.output_name = 'report'
            print("Output name:   " + str(args.output_name))

        if not os.path.exists(args.output):
            os.makedirs(args.output)


        self.imagePath = makePath(args.output) + args.output_name
        if not os.path.exists(self.imagePath):
            os.makedirs(self.imagePath)

        pltcfg = PlotConfig()
        pltcfg.saveToFile(self.imagePath)

        args.pltcfg = pltcfg

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)


    def execParallel(self, data, environment):
        iFilesInFolder = 0

        localEnv = {}

        for report in self.dReporters:
            reportObj = self.dReporters[report]
            localEnv[report] = reportObj._createLocalEnvironment()


        f5folder = Fast5Directory(data)

        for file in f5folder.collect():

            for report in self.dReporters:

                reportObj = self.dReporters[report]
                reportObj.handleEntity(file, localEnv[report], environment)

            iFilesInFolder += 1


        print("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return localEnv

    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        for report in self.dReporters:

            if not report in existResult:
                existResult[report] = None

            reportObj = self.dReporters[report]
            existResult[report] = reportObj.joinParallel( existResult[report], newResult[report], oEnvironment )

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):


        with open(args.output + args.output_name + ".html", 'w') as htmlFile:

            htmlFile.write("<html><body>\n")

            for report in self.dReporters:

                htmlFile.write("<h1>" + str(report) + "</h1>\n")

                reportObj = self.dReporters[report]

                reportObj.makeResults(parallelResult[report], oEnvironment, args)

                createdPlots = args.pltcfg.getCreatedPlots()
                args.pltcfg.resetCreatedPlots()

                for x in createdPlots:

                    relPath = os.path.relpath( x, args.output )

                    htmlFile.write("<p><img src=\"" + str(relPath) + "\"/></p>\n")
                    print(str(report) + "\t" + str(relPath))

            htmlFile.write("</body></html>\n")



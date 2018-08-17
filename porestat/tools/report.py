import argparse
import os

from porestat.tools.experiment_ls import ExperimentLs
from porestat.tools.kmer_coverage import KmerHistogram
from .quality_position import QualityPosition
from .length_histogram import LengthHistogram
from .nucleotide_distribution import NucleotideDistribution
from .quality_distribution import QualityDistribution
from ..utils.Files import makePath
from .channel_occupancy import ChannelOccupancy
from .yield_plot import YieldPlot


from ..plots.plotconfig import PlotConfig, PlotSaveTYPE
from ..hdf5tool.Fast5File import Fast5TYPEAction

from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

from collections import OrderedDict
import dill as pickle

class ReportFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(ReportFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser = subparsers.add_parser('summary', help='read performance summary')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)

        parser.add_argument('-o', '--output', type=str, help='output folder for report', required=True)
        parser.add_argument('-n', '--output-name', type=str, help='output name', required=False)

        parser.add_argument('--no-read-type-subplot', dest='addTypeSubplot', action='store_false', default=True, help='do not add type subplots')
        parser.add_argument('-q', '--read-type', nargs='+', dest='read_type', action=Fast5TYPEAction, default=None)
        parser.add_argument('-u', '--user-run', dest='user_run', action='store_true', default=False)

        parser.add_argument('--reference', dest='reference', type=argparse.FileType('r'), help='fasta reference for kmer distribution', default=None)


        parser.add_argument('--save-parallel-result', type=str, default=None)
        parser.add_argument('--load-parallel-result', nargs='+', type=str, default=None, help='specify any saved pickle files to combine for report')


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
            ('OVERVIEW', ExperimentLs(args)),
            ('OCCUPANCY', ChannelOccupancy(args)),
            ('YIELD', YieldPlot(args)),
            ('LENGTH', LengthHistogram(args)),
            ('NUCLEOTIDES', NucleotideDistribution(args)),
            ('QUALITY DISTRIBUTION', QualityDistribution(args)),
            ('QUALITY BY POSITION', QualityPosition(args)),
            ('KMER DISTRIBUTION', KmerHistogram(args))
        ])

        #self.dReporters = OrderedDict([('OVERVIEW', ExperimentLs(args)), ('QUALITY BY POSITION', QualityPosition(args))])

        self.dReportersArgs = OrderedDict([

            ('YIELD', {
                'separate_subplots': True
            }),
            ('KMER DISTRIBUTION', {
                'violin': False,
                'k': 5,
                'mc': 10
            })

        ])

        args.output = makePath(args.output)

        print("Output folder: " + str(args.output))
        print("Output name:   " + str(args.output_name))

        if args.output_name == None:
            args.output_name = 'report'
            print("Output name:   " + str(args.output_name))

        if not os.path.exists(args.output):
            os.makedirs(args.output)

        self.data_path = makePath(makePath(args.output) + args.output_name)
        if not os.path.exists(self.data_path):
            os.makedirs(self.data_path)

        import mpld3, shutil
        d3js_path = mpld3.getD3js()
        mpld3_path = mpld3.getmpld3js(True)

        d3js_dest = self.data_path + os.path.split(d3js_path)[1]
        mpld3_dest = self.data_path + os.path.split(mpld3_path)[1]

        shutil.copyfile(d3js_path, d3js_dest)
        shutil.copyfile(mpld3_path, mpld3_dest)

        self.html_path = args.output + args.output_name + ".html"

        args.pltcfg.d3js = os.path.relpath(d3js_dest, args.output)
        args.pltcfg.mpld3js = os.path.relpath(mpld3_dest, args.output)

        print("Relative js path: " + args.pltcfg.d3js)


    def prepareInputs(self, args):

        if args.load_parallel_result == None:
            return self.manage_folders_reads(args)
        else:
            llResultExists = all([os.path.isfile(x) for x in args.load_parallel_result])

            if llResultExists:
                return []

        return self.manage_folders_reads(args)



    def execParallel(self, data, environment):
        iFilesInFolder = 0

        localEnv = {}

        for folder in data:

            for report in self.dReporters:
                reportObj = self.dReporters[report]
                localEnv[report] = reportObj._createLocalEnvironment()

            f5folder = Fast5Directory(folder)

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

        args.pltcfg.setOutputType(PlotSaveTYPE.HTML_STRING)

        if not args.save_parallel_result is None:

            with open(args.save_parallel_result, 'wb') as pickleFile:
                pickle.dump( parallelResult, pickleFile )

            print("Result saved to: " + args.save_parallel_result)

        if not args.load_parallel_result is None and all([os.path.isfile(x) for x in args.load_parallel_result]):
            print("Trying to load result from: " + ", ".join(args.load_parallel_result))

            for res in args.load_parallel_result:
                with open(res, 'rb') as pickleFile:
                    loadedResult = pickle.load(pickleFile)

                    if parallelResult == None:
                        parallelResult = loadedResult
                    else:
                        parallelResult = mergeDicts(parallelResult, loadedResult)

                print("Result loaded: " + res)

        with open(args.output + args.output_name + ".html", 'w') as htmlFile:

            mpld3js = "<script src=" + args.pltcfg.mpld3js +"></script>\n"
            d3js = "<script src=" + args.pltcfg.d3js + "></script>\n"

            htmlFile.write(
                """
                <html>
                <head>
                """
                +d3js+
                mpld3js+
                """
                
                <style>
                * {
                    font-family: "Arial", Verdana, sans-serif;
                    }
                </style>
                </head>
                
                
                </head>
                <body>
                """
            )

            for report in self.dReporters:

                print("Running report: " + str(report))

                reporterArgs = self.prepareEnvironment(args)
                reporterArgs.output = None
                reporterArgs.output_type = None
                reporterArgs.pltcfg = args.pltcfg

                reporterArgs = self.patchArgs(reporterArgs, report)

                reporterArgs.pltcfg.saveToFile(self.data_path + "/" + report)

                htmlFile.write("<h1>" + str(report) + "</h1>\n")

                reportObj = self.dReporters[report]

                reportObj.makeResults(parallelResult[report], oEnvironment, reporterArgs)
                args.pltcfg = reporterArgs.pltcfg

                createdPlots = args.pltcfg.getCreatedPlots()
                args.pltcfg.resetCreatedPlots()

                for x in createdPlots:

                    relPath = os.path.relpath( x, args.output )
                    htmlFile.write(x + "\n")

                    #print(str(report) + "\t" + str(relPath))

                htmlFile.flush()

            htmlFile.write("</body></html>\n")


    def patchArgs(self, args, reporter):

        if not reporter in self.dReportersArgs:
            return args

        repArgs =self.dReportersArgs[reporter]
        if repArgs == None or len(repArgs) == 0:
            return args

        for x in repArgs:
            args.__dict__[x] = repArgs[x]

        return args
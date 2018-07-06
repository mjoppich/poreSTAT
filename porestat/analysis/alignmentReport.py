import argparse
import os

import HTSeq
from porestat.utils.ArgParseExt import FileStubType, FolderType

from .alignmentStatistics import AlignmentStatisticAnalysis
from .read_counts import ReadCountAnalysis

from ..utils.Files import makePath

from ..plots.plotconfig import PlotConfig, PlotSaveTYPE
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter, defaultdict
from ..utils.Utils import mergeDicts, mergeCounter

from collections import OrderedDict
import dill as pickle




class ReportFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(ReportFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-s', '--sam', nargs='+', type=argparse.FileType('r'), help='alignment files', required=True)
        parser.add_argument('-f', '--fasta', dest='fasta_files', nargs='+', type=argparse.FileType('r'), help='read inputs for alignment', required=True)
        parser.add_argument('-r', '--read-info', nargs='+', type=argparse.FileType('r'), help='read summary file', required=False)

        parser.add_argument('-o', '--output', type=FolderType('w'), help='output folder for report', required=True)
        parser.add_argument('-n', '--output-name', type=str, help='output name', required=False)

        parser.add_argument('--save-parallel-result', type=argparse.FileType('w'), default=None)
        parser.add_argument('--load-parallel-result', nargs='+', type=str, default=None, help='specify any saved pickle files to combine for report')

        parser.add_argument('-v', '--violin', dest='violin', action='store_true', default=False)


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

        self.args.pltcfg.setOutputType(PlotSaveTYPE.HTML_STRING)



        self.dReporters = OrderedDict([
            ('ALIGNMENT', AlignmentStatisticAnalysis(args)),
            #('READ COUNTS/COVERAGE', ReadCountAnalysis(args))
        ])

        self.dReportersArgs = OrderedDict([
            ('ALIGNMENT', {
                'violin': args.violin,
                'mc': 50,
                'errork': 5,
                'perfectk': 21
            })
        ])

        self.filesPerReport = {}
        self.reporterEnvironments = {}

    def prepareInputs(self, args):

        self.filesPerReport = {}
        for report in self.dReporters:
            reportObj = self.dReporters[report]
            self.filesPerReport[report] = reportObj.prepareInputs(args)

        allFiles = set()
        for reporter in self.filesPerReport:
            for file in self.filesPerReport[reporter]:
                allFiles.add(file)

        if args.load_parallel_result != None:
            llResultExists = all([os.path.isfile(x) for x in args.load_parallel_result])

            if llResultExists:
                return []

        return list(allFiles)

    def prepareEnvironment(self, args):

        self.reporterEnvironments = {}
        for report in self.dReporters:
            reportObj = self.dReporters[report]
            self.reporterEnvironments[report] = reportObj.prepareEnvironment(args)

        return self._makeArguments(args)


    def execParallel(self, data, environment):

        retTuples = []
        for alignedFile in data:

            print("Processing File:", alignedFile.name)

            if alignedFile.name.endswith(".bam"):
                opener = HTSeq.BAM_Reader
            else:
                opener = HTSeq.SAM_Reader

            iProcessedAlignments = 0
            localEnv = {}
            for report in self.dReporters:
                reportObj = self.dReporters[report]
                localEnv[report] = reportObj._createLocalEnvironment()


            for readAlignment in opener(alignedFile):

                for report in self.dReporters:
                    if alignedFile in self.filesPerReport[report]:
                        reportObj = self.dReporters[report]
                        localEnv[report] = reportObj.handleEntity(readAlignment, localEnv[report], self.reporterEnvironments[report])

                iProcessedAlignments += 1

            print("File " + str(alignedFile.name) + ": " + str(iProcessedAlignments) + "alignments processed")

            retTuples.append( (alignedFile.name, localEnv) ) #(data, localEnv)

        return retTuples


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        for elem in newResult:
            (file, newResults) = elem

            for report in newResults:

                if not report in existResult:
                    existResult[report] = None

                reportObj = self.dReporters[report]
                existResult[report] = reportObj.joinParallel( existResult[report], [(file, newResults[report])], self.reporterEnvironments[report] )

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):


        args.output = makePath(args.output)

        print("Output folder: " + str(args.output))
        print("Output name:   " + str(args.output_name))

        if args.output_name == None:
            args.output_name = 'report'
            print("Changed Output name:   " + str(args.output_name))

        for report in self.dReporters:
            print("Running report: " + str(report))

            reporterArgs = self.prepareEnvironment(args)
            reporterArgs.output = None
            reporterArgs.output_type = None
            reporterArgs.pltcfg = args.pltcfg

            reporterArgs = self.patchArgs(reporterArgs, report)

            reporterArgs.pltcfg.saveToFile(args.output + "/" + report)

            reporterArgs.pltcfg.addHTMLPlot("<h1>" + str(report) + "</h1>\n")

            reportObj = self.dReporters[report]

            reportObj.makeResults(parallelResult[report], oEnvironment, reporterArgs)
            args.pltcfg = reporterArgs.pltcfg

        args.pltcfg.prepareHTMLOutput(args.output, args.output_name + ".html", relativeImport=True)


    def patchArgs(self, args, reporter):

        if not reporter in self.dReportersArgs:
            return args

        repArgs =self.dReportersArgs[reporter]
        if repArgs == None or len(repArgs) == 0:
            return args

        for x in repArgs:
            args.__dict__[x] = repArgs[x]

        return args
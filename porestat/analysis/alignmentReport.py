import argparse
import os
import time

import HTSeq
from porestat.analysis.ParallelAlignmentPSTReportableInterface import ParallelAlignmentPSTReportableInterface

from porestat.utils.Parallel import MapReduce

from porestat.utils.ArgParseExt import FileStubType, FolderType

from .alignmentStatistics import AlignmentStatisticAnalysis
from .read_counts import ReadCountAnalysis

from ..utils.Files import makePath, eprint

from ..plots.plotconfig import PlotConfig, PlotSaveTYPE
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter, defaultdict
from ..utils.Utils import mergeDicts, mergeCounter

from collections import OrderedDict
import dill as pickle


import pysam

#python3 scripts/poreAnalysis.py report --load-parallel-result /mnt/g/sequ_into_demo/9_COVID_RNA/pass.3.pickle --sam /mnt/g/sequ_into_demo/9_COVID_RNA/reads.sam --read-info /mnt/g/sequ_into_demo/9_COVID_RNA/reads.info --gtf /mnt/g/sequ_into_demo/9_COVID_RNA/Sars_cov_2.ASM985889v3.100.gtf --fasta /mnt/g/sequ_into_demo/9_COVID_RNA/Sars_cov_2.ASM985889v3.dna.toplevel.fa --output /mnt/g/sequ_into_demo/9_COVID_RNA/sam_analysis

class ReportFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(ReportFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-s', '--sam', nargs='+', type=argparse.FileType('r'), help='alignment files', required=True)
        parser.add_argument('-f', '--fasta', dest='fasta_files', nargs='+', type=argparse.FileType('r'), help='read inputs for alignment', required=True)
        parser.add_argument('-r', '--read-info', nargs='+', type=argparse.FileType('r'), help='read summary file', required=False)
        parser.add_argument('-g', '--gtf', type=argparse.FileType('r'), help="gtf file for features",
                            required=False, default=None)

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


class ReportAnalysis(ParallelAlignmentPSTReportableInterface):

    def __init__(self, args):

        super(ReportAnalysis, self).__init__( args )

        self.args.pltcfg.setOutputType(PlotSaveTYPE.HTML_STRING)
        self.chunkSize = 10


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

        localEnv = {}
        for report in self.dReporters:
            reportObj = self.dReporters[report]
            localEnv[report] = reportObj._createLocalEnvironment()

        return self._makeArguments(args), localEnv


    def handleEntity(self, samHeader, stringReads, envs ):

        alignmentFileName = envs[0]
        globalEnv = envs[1]
        localEnv = envs[2]


        for samReadAlignment in stringReads:

            aln = pysam.AlignedSegment.fromstring(str, samReadAlignment, samHeader)

            for report in self.dReporters:
                if alignmentFileName in self.filesPerReport[report]:
                    reportObj = self.dReporters[report]
                    localEnv[report] = reportObj.handleEntity(aln, localEnv[report], self.reporterEnvironments[report])

        return [(alignmentFileName, localEnv)]


    def joinParallel(self, existResult, newResult, oEnvironment):

        print("Report Join Parallel Starts")


        if existResult == None:
            existResult = {}

        for elem in newResult:
            (file, newResults) = elem

            print(len(newResults))

            for report in newResults:

                if not report in existResult:
                    existResult[report] = None

                reportObj = self.dReporters[report]
                existResult[report] = reportObj.joinParallel( existResult[report], [(file, newResults[report])], self.reporterEnvironments[report] )

        print("Report Join Parallel Ends")

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

            reporterArgs, localEnv = self.prepareEnvironment(args)
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
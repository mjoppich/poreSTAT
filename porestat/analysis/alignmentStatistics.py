import argparse

import pickle
from collections import defaultdict

from .ParallelAlignmentPSTReportableInterface import ParallelAlignmentPSTReportableInterface
from ..hdf5tool.Fast5File import Fast5TYPEAction
from ..plots.plotconfig import PlotConfig
from ..hdf5tool import Fast5TYPE
from ..plots.poreplot import PorePlot

from ..utils.Files import fileExists
from ..utils.DataFrame import DataFrame

from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException
from ..utils.Utils import mergeDicts
from ..utils.Stats import calcN50

from collections import Counter
import HTSeq

class MinimalSeq:

    def __init__(self, seq, name, desc):

        self.seq = seq
        self.name = name
        self.desc = desc

        self.htseq = None

    def toHTSeq(self):

        if self.htseq == None:

            htseq = HTSeq.Sequence(self.seq, self.name)
            htseq.descr = self.desc

            self.htseq = htseq

        return self.htseq

class AlignmentStatisticAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(AlignmentStatisticAnalysisFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('align_stats', help='Alignment statistics (amount aligned, indels, ...)')
        parser.add_argument('-s', '--sam', nargs='+', type=str, help='alignment files', required=True)
        parser.add_argument('-f', '--fasta', dest='fasta_files', nargs='+', type=str, help='read inputs for alignment', required=True)
        parser.add_argument('-r', '--read-info', nargs='+', type=str, help='read summary file', required=False)

        parser.add_argument('-q', '--read-type', nargs='+', dest='read_type', action=Fast5TYPEAction, default=None)
        parser.add_argument('-v', '--violin', dest='violin', action='store_true', default=False)

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return AlignmentStatisticAnalysis(simArgs)



class AlignmentStatisticAnalysis(ParallelAlignmentPSTReportableInterface):

    def __init__(self, args):

        super(AlignmentStatisticAnalysis, self).__init__(args)

        self.cigar_ops = ['M', 'X', '=', 'D', 'I', 'S', 'X', 'H', 'N', 'P']

        self.readInfoIdx = {}


    def readReadInfo(self, args):

        for readInfoFile in args.read_info:

            if not fileExists(readInfoFile):
                raise PSToolException("Read info file does not exist: " + str(readInfoFile))

        args.readInfo = None

        for readInfoFile in args.read_info:
            allReadData = DataFrame.parseFromFile(readInfoFile, cDelim='\t')

            samIdx = allReadData.addColumn('SAM_READ_NAME')
            readNameIdx = allReadData.getColumnIndex('READ_NAME')
            allReadData.applyByRow(samIdx, lambda x: x[readNameIdx].split(' ')[0])

            if args.readInfo == None:
                args.readInfo = allReadData
            else:
                args.readInfo.merge( allReadData )

    def readSequences(self, args):

        for fastaFile in args.fasta_files:

            if not fileExists(fastaFile):
                raise PSToolException("Fasta file does not exist: " + str(fastaFile))


        if not self.hasArgument('fasta', args) or args.fasta == None:
            args.fasta = {}

            for fastaFile in args.fasta_files:

                for seq in HTSeq.FastaReader( fastaFile ):
                    args.fasta[seq.name] = MinimalSeq(seq.seq, seq.name, seq.descr)

    def _makePropDict(self):

        props = defaultdict(list)
        return props

    def prepareInputs(self, args):

        self.readReadInfo(args)
        self.readSequences(args)

        for x in args.sam:
            if not fileExists(x):
                PSToolException("sam file does not exist: " + str(x))

        return [x for x in args.sam]

    def _createLocalEnvironment(self):
        return { 'BASE_COMPOSITION': defaultdict(Counter),
                 'UNALIGNED_READS': set(),
                 'READ_STATS': {}
                 }

    def handleEntity(self, readAlignment, localEnv, globalEnv):

        readName = readAlignment.read.name
        readInfo = globalEnv.readInfo.findRow('SAM_READ_NAME', readName)
        if not readInfo == None:
            readType = Fast5TYPE[ readInfo[ "TYPE" ] ]
            readID = readInfo[ "READ_ID" ]
        else:
            print("read id none for " + str(readName))

        if readAlignment.aligned:

            readInterval = readAlignment.iv
            cigarOfRead = readAlignment.cigar

            cigarOp2Length = self._makePropDict()

            for cigarOp in cigarOfRead:

                readSeq = readAlignment.read_as_aligned[cigarOp.query_from:cigarOp.query_to]

                for base in readSeq.seq:
                    localEnv['BASE_COMPOSITION'][cigarOp.type][chr(base)] += 1

                if cigarOp.type == 'M':

                    fastaReq = globalEnv.fasta[ readInterval.chrom ]

                    fastaReq = fastaReq.toHTSeq()

                    fastaSeq = fastaReq[ cigarOp.ref_iv.start:cigarOp.ref_iv.end ]

                    editDistance = self.calculateEditDistance( readSeq, fastaSeq )

                    cigarOp2Length['M'].append(cigarOp.size)
                    cigarOp2Length['X'].append(editDistance)
                    cigarOp2Length['='].append(cigarOp.size - editDistance)

                else:
                    cigarOp2Length[ cigarOp.type ].append(cigarOp.size)

            localEnv['READ_STATS'][(readID, readType)] = cigarOp2Length

        else:
            localEnv['UNALIGNED_READS'].add(readID)

        return localEnv

    def calculateEditDistance(self, htSeq1, htSeq2):

        seq1 = htSeq1.seq
        seq2 = htSeq2.seq

        diff = 0
        for i in range(0, min(len(seq1), len(seq2))):
            if seq1[i] != seq2[i]:
                diff += 1

        return diff


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        (file, data) = newResult

        if file in existResult:
            existResult[file] = mergeDicts(existResult[file], data)
        else:
            existResult[file] = data

        return existResult


    def makeResults(self, allFilesResult, oEnvironment, args):

        for fileName in allFilesResult:

            parallelResult = allFilesResult[fileName]

            readCIGARs = parallelResult['READ_STATS']
            alignedIDs = set([x[0] for x in readCIGARs])

            unalignedIDs = parallelResult['UNALIGNED_READS']
            additionalUnalignedReads = set()

            for row in args.readInfo:

                readid = row['READ_ID']

                if readid in alignedIDs:
                    continue

                if readid in unalignedIDs:
                    continue

                additionalUnalignedReads.add(readid)

            plotDataDict = {}
            plotDataDict['Aligned'] = len(alignedIDs)
            plotDataDict['Unaligned'] = len(unalignedIDs)
            plotDataDict['Not Aligned'] = len(additionalUnalignedReads)

            plotDict = {'All Reads': plotDataDict}

            PorePlot.plotBars(plotDict, fileName + ": Alignment Result", "", "Count", pltcfg=args.pltcfg)

            overviewByAligned = defaultdict(Counter)
            for row in args.readInfo:

                readid = row['READ_ID']
                read_type = row['TYPE']

                if readid in alignedIDs:
                    overviewByAligned['ALIGNED'][read_type] += 1
                elif readid in unalignedIDs:
                    overviewByAligned['UNALIGNED'][read_type] += 1
                elif readid in additionalUnalignedReads:
                    overviewByAligned['NOT ALIGNED'][read_type] += 1
                else:
                    overviewByAligned['UNKNOWN'][read_type] += 1

            print(overviewByAligned)

            PorePlot.plotBars(overviewByAligned, fileName + ": Alignment Result By Type", "", "Count", pltcfg=args.pltcfg)

            if not self.hasArgument('read_type', args) or not args.read_type:

                totalCounter = defaultdict(list)

                for counterPair in readCIGARs:

                    readID = counterPair[0]
                    readType = counterPair[1]
                    readCounter = readCIGARs[counterPair]

                    for x in readCounter:
                        totalCounter[x] += readCounter[x]

                if self.hasArgument('violin', args) and args.violin:
                    PorePlot.plotViolin(totalCounter, [x for x in totalCounter], "CIGARs: all types", pltcfg=args.pltcfg)
                else:
                    PorePlot.plotBoxplot(totalCounter, [x for x in totalCounter], "CIGARs: all types", pltcfg=args.pltcfg)

            totalCounter = defaultdict(lambda: defaultdict(list))

            for counterPair in readCIGARs:

                readID = counterPair[0]
                readType = counterPair[1]
                readCounter = readCIGARs[counterPair]

                for x in readCounter:
                    totalCounter[readType][x] += readCounter[x]

            for readType in totalCounter:
                if self.hasArgument('violin', args) and args.violin:
                    PorePlot.plotViolin(totalCounter[readType], [x for x in totalCounter[readType]], "CIGARs: " + str(readType),
                                        pltcfg=args.pltcfg)
                else:
                    PorePlot.plotBoxplot(totalCounter[readType], [x for x in totalCounter[readType]], "CIGARs: " + str(readType),
                                         pltcfg=args.pltcfg)
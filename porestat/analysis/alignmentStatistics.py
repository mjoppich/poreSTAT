import argparse

import pickle
from collections import defaultdict

from porestat.hdf5tool.Fast5File import Fast5TYPEAction
from ..plots.plotconfig import PlotConfig
from porestat.hdf5tool import Fast5TYPE
from porestat.plots.poreplot import PorePlot

from ..utils.Files import fileExists
from ..utils.DataFrame import DataFrame
from ..tools.ParallelPTTInterface import ParallelPSTInterface
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
        parser.add_argument('-f', '--fasta', nargs='+', type=str, help='read inputs for alignment', required=True)
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



class AlignmentStatisticAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(AlignmentStatisticAnalysis, self).__init__(args)

        self.readInfo = None
        self.fasta = None
        self.genome_annotation = None

        self.cigar_ops = ['M', 'X', '=', 'D', 'I', 'S', 'X', 'H', 'N', 'P']

        self.readInfoIdx = {}


    def readReadInfo(self, args):

        for readInfoFile in args.read_info:

            if not fileExists(readInfoFile):
                raise PSToolException("Read info file does not exist: " + str(readInfoFile))

        self.readInfo = None

        for readInfoFile in args.read_info:
            allReadData = DataFrame.parseFromFile(readInfoFile, cDelim='\t')

            samIdx = allReadData.addColumn('SAM_READ_NAME')
            readNameIdx = allReadData.getColumnIndex('READ_NAME')
            allReadData.applyByRow(samIdx, lambda x: x[readNameIdx].split(' ')[0])

            if self.readInfo == None:
                self.readInfo = allReadData
            else:
                self.readInfo.merge( allReadData )

    def readSequences(self, args):

        for fastaFile in args.fasta:

            if not fileExists(fastaFile):
                raise PSToolException("Fasta file does not exist: " + str(fastaFile))

        self.fasta = {}

        for fastaFile in args.fasta:

            for seq in HTSeq.FastaReader( fastaFile ):
                self.fasta[seq.name] = MinimalSeq(seq.seq, seq.name, seq.descr)

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

    def prepareEnvironment(self, args):

        env = self._makeArguments(args)
        env.fasta = self.fasta

        return env

    def execParallel(self, data, env):

        if data.endswith(".bam"):
            opener = HTSeq.BAM_Reader
        else:
            opener = HTSeq.SAM_Reader


        baseCompositions = defaultdict(Counter)
        allReadStats = {}

        unalignedIDs = set()

        iProcessedReads = 0

        for readAlignment in opener( data ):

            readName = readAlignment.read.name
            readType = Fast5TYPE.UNKNOWN

            iProcessedReads += 1

            if iProcessedReads % 1000 == 0:
                print("Done: " + str(iProcessedReads))

            readInfo = self.readInfo.findRow('SAM_READ_NAME', readName)
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
                        baseCompositions[cigarOp.type][chr(base)] += 1

                    if cigarOp.type == 'M':

                        fastaReq = env.fasta[ readInterval.chrom ]

                        fastaReq = fastaReq.toHTSeq()

                        fastaSeq = fastaReq[ cigarOp.ref_iv.start:cigarOp.ref_iv.end ]

                        editDistance = self.calculateEditDistance( readSeq, fastaSeq )

                        cigarOp2Length['M'].append(cigarOp.size)
                        cigarOp2Length['X'].append(editDistance)
                        cigarOp2Length['='].append(cigarOp.size - editDistance)

                    else:
                        cigarOp2Length[ cigarOp.type ].append(cigarOp.size)

                allReadStats[readID] = cigarOp2Length

            else:
                unalignedIDs.add(readID)





        return ( allReadStats, unalignedIDs )

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
            existResult = ({}, set())

        returnResult = (mergeDicts(existResult[0], newResult[0]), existResult[1].union(newResult[1]))
        return returnResult


    def makeResults(self, parallelResult, oEnvironment, args):

        readCIGARs = parallelResult[0]
        alignedIDs = set([x for x in parallelResult[0]])

        unalignedIDs = parallelResult[1]

        additionalUnalignedReads = set()

        for row in self.readInfo:

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

        print(plotDict)

        PorePlot.plotBars(plotDict, "Alignment Stat", "X", "Y", pltcfg=args.pltcfg)


        overviewByAligned = defaultdict(Counter)
        for row in self.readInfo:

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

        PorePlot.plotBars(overviewByAligned, "Alignment Stat", "X", "Y", pltcfg=args.pltcfg)

        if not args.read_type:

            totalCounter = defaultdict(list)

            for readID in readCIGARs:

                cigarCounter = readCIGARs[readID]

                for x in cigarCounter:
                    totalCounter[x] += cigarCounter[x]

            if args.violin:
                PorePlot.plotViolin(totalCounter, [x for x in totalCounter], "CIGARs", pltcfg=args.pltcfg)
            else:
                PorePlot.plotBoxplot(totalCounter, [x for x in totalCounter], "CIGARs", pltcfg=args.pltcfg)

        else:
            totalCounter = defaultdict( defaultdict(list) )

            for counterPair in readCIGARs:

                readType = counterPair[0]
                readCounter = counterPair[1]

                for x in readCounter:
                    totalCounter[readType][x] += readCounter

            for readType in totalCounter:
                if args.violin:
                    PorePlot.plotViolin(totalCounter[readType], [x for x in totalCounter[readType]], "CIGARs", pltcfg=args.pltcfg)
                else:
                    PorePlot.plotBoxplot(totalCounter[readType], [x for x in totalCounter[readType]], "CIGARs", pltcfg=args.pltcfg)








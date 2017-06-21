import argparse

import pickle
from collections import defaultdict

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

class BasecallFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(BasecallFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('align_stats', help='Alignment statistics (amount aligned, indels, ...)')
        parser.add_argument('-s', '--sam', nargs='+', type=str, help='alignment files', required=True)
        parser.add_argument('-f', '--fasta', nargs='+', type=str, help='read inputs for alignment', required=True)
        parser.add_argument('-r', '--read-info', nargs='+', type=str, help='read summary file', required=False)

        parser.add_argument('-q', '--read-type', nargs='+', type=str, choices=[x for x in Fast5TYPE.str2type], help='read types ('+ ",".join([x for x in Fast5TYPE.str2type]) +')')
        parser.add_argument('-v', '--violin', dest='violin', action='store_true', default=False)

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return BasecallFactory(simArgs)



class Basecall(ParallelPSTInterface):

    def __init__(self, args):

        super(Basecall, self).__init__(args)

        self.readInfo = None
        self.fasta = None
        self.genome_annotation = None

        self.cigar_ops = ['M', 'X', '=', 'D', 'I', 'S', 'X', 'H', 'N', 'P']


    def readReadInfo(self, args):

        for readInfoFile in args.read_info:

            if not fileExists(readInfoFile):
                raise PSToolException("Read info file does not exist: " + str(readInfoFile))

        self.readInfo = {}

        for readInfoFile in args.read_info:
            allReadData = DataFrame.parseFromFile(readInfoFile, cDelim='\t')

            readNameIdx = allReadData.getColumnIndex("READ_NAME")
            readTypeIdx = allReadData.getColumnIndex("TYPE")
            readRunIdx = allReadData.getColumnIndex("USER_RUN_NAME")

            for elem in allReadData.vElements:

                readName = elem[readNameIdx].split(" ")[0]
                readType = elem[readTypeIdx]
                readRun  = elem[readRunIdx]

                self.readInfo[readName] = (readType, readRun)

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

        allReadStats = []

        for readAlignment in opener( data ):

            if readAlignment.aligned:

                readName = readAlignment.read.name

                readType = Fast5TYPE.UNKNOWN
                if readName in self.readInfo:
                    readType = self.readInfo[readName][0]

                readInterval = readAlignment.iv
                cigarOfRead = readAlignment.cigar

                cigarOp2Length = self._makePropDict()

                for cigarOp in cigarOfRead:

                    if cigarOp.type == 'M':

                        fastaReq = env.fasta[ readInterval.chrom ]

                        fastaReq = fastaReq.toHTSeq()

                        fastaSeq = fastaReq[ cigarOp.ref_iv.start:cigarOp.ref_iv.end ]
                        readSeq = readAlignment.read_as_aligned[cigarOp.query_from:cigarOp.query_to]

                        editDistance = self.calculateEditDistance( readSeq, fastaSeq )

                        cigarOp2Length['M'].append(cigarOp.size)
                        cigarOp2Length['X'].append(editDistance)
                        cigarOp2Length['='].append(cigarOp.size - editDistance)

                    else:
                        cigarOp2Length[ cigarOp.type ].append(cigarOp.size)

                print("done read: " + readName)
                allReadStats.append( (readType, cigarOp2Length) )

        return allReadStats

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
            existResult = []

        existResult += newResult

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):


        if not args.read_type:

            totalCounter = defaultdict(list)

            for counterPair in parallelResult:

                for x in counterPair[1]:
                    totalCounter[x] += counterPair[1][x]

            if args.violin:
                PorePlot.plotViolin(totalCounter, [x for x in totalCounter], "CIGARs", pltcfg=args.pltcfg)
            else:
                PorePlot.plotBoxplot(totalCounter, [x for x in totalCounter], "CIGARs", pltcfg=args.pltcfg)

        else:
            totalCounter = defaultdict( defaultdict(list) )

            for counterPair in parallelResult:

                readType = counterPair[0]
                readCounter = counterPair[1]

                for x in readCounter:
                    totalCounter[readType][x] += readCounter

            for readType in totalCounter:
                if args.violin:
                    PorePlot.plotViolin(totalCounter[readType], [x for x in totalCounter[readType]], "CIGARs", pltcfg=args.pltcfg)
                else:
                    PorePlot.plotBoxplot(totalCounter[readType], [x for x in totalCounter[readType]], "CIGARs", pltcfg=args.pltcfg)








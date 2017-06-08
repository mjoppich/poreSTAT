import argparse

import pickle

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

class AlignmentStatisticAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(AlignmentStatisticAnalysisFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('align_stats', help='Alignment statistics (amount aligned, indels, ...)')
        parser.add_argument('-s', '--sam', nargs='+', type=str, help='alignment files', required=True)
        parser.add_argument('-f', '--fasta', nargs='+', type=str, help='read inputs for alignment', required=True)
        parser.add_argument('-r', '--read-info', nargs='+', type=str, help='read summary file', required=False)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return AlignmentStatisticAnalysis(simArgs)



class AlignmentStatisticAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(AlignmentStatisticAnalysis, self).__init__(args)

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

                readName = elem[readNameIdx]
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

        props = {}

        for cigarOp in self.cigar_ops:
            props[cigarOp] = []

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

        copy = pickle.dump(self.fasta, open("/tmp/bla", 'wb'))

        return env

    def execParallel(self, data, env):

        if data.endswith(".bam"):
            opener = HTSeq.BAM_Reader
        else:
            opener = HTSeq.SAM_Reader

        for readAlignment in opener( data ):

            if readAlignment.aligned:

                readName = readAlignment.read.name

                readInterval = readAlignment.iv
                cigarOfRead = readAlignment.cigar

                cigarOp2Length = self._makePropDict()

                for cigarOp in cigarOfRead:

                    if cigarOp.type == 'M':

                        fastaReq = env.fasta[ readInterval.chrom ]

                        fastaReq = HTSeq.Sequence(fastaReq.seq, fastaReq.name)
                        fastaSeq = fastaReq[ cigarOp.ref_iv.start:cigarOp.ref_iv.end ]

                        if readInterval.strand == '-':
                            fastaSeq = fastaSeq.get_reverse_complement()

                        editDistance = self.calculateEditDistance( readAlignment.read[cigarOp.query_from:cigarOp.query_to], fastaSeq )

                        cigarOp2Length['M'].append(cigarOp.size)
                        cigarOp2Length['X'].append(editDistance)
                        cigarOp2Length['='].append(cigarOp.size - editDistance)

                    else:
                        cigarOp2Length[ cigarOp.type ].append(cigarOp.size)

        return None

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

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'FILES', 'AVG_LENGTH', 'N50', 'BASECALL_2D', 'BASECALL_1D_COMPL',
                            'BASECALL_1D', 'BASECALL_RNN_1D', 'PRE_BASECALL', 'UNKNOWN']

        if parallelResult == None:
            raise PSToolException("No Files collected")

        allobservations = {}
        for runid in parallelResult:

            allLengths = []
            countByType = Counter()

            for type in parallelResult[runid]['TYPE']:

                lengths = parallelResult[runid]['TYPE'][type]

                allLengths += lengths

                countByType[type] += len(lengths)

            fileCount = len(allLengths)
            avgLength = sum(allLengths) / fileCount
            (n50, l50) = calcN50(allLengths)
            run_user_name = parallelResult[runid]['NAMES']['USER_RUN_NAME']

            observations = {

                'RUNID': runid,
                'USER_RUN_NAME': ",".join(run_user_name),
                'FILES': fileCount,
                'AVG_LENGTH': avgLength,
                'N50': n50,
                'BASECALL_2D': countByType[self.str2fileType['BASECALL_2D']],
                'BASECALL_1D_COMPL': countByType[self.str2fileType['BASECALL_1D_COMPL']],
                'BASECALL_1D': countByType[self.str2fileType['BASECALL_1D']],
                'BASECALL_RNN_1D': countByType[self.str2fileType['BASECALL_RNN_1D']],
                'PRE_BASECALL': countByType[self.str2fileType['PRE_BASECALL']],
                'UNKNOWN': countByType[self.str2fileType['UNKNOWN']]
            }

            allobservations[runid] = observations


        sortedruns = sorted([x for x in allobservations])

        print("\t".join(makeObservations))

        for runid in sortedruns:

            allobs = []
            for x in makeObservations:
                allobs.append(str(allobservations[runid][x]))

            print("\t".join(allobs))



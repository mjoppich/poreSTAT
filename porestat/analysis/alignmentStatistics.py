import argparse

from ..utils.DataFrame import DataFrame
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Parallel import Parallel as ll
from ..utils.Utils import mergeDicts, fileExists
from ..utils.Stats import calcN50

import HTSeq

class AlignmentStatisticAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(AlignmentStatisticAnalysisFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('align_stats', help='Alignment statistics (amount aligned, indels, ...)')
        parser.add_argument('-s', '--sam', nargs='+', type=argparse.FileType('r'), help='alignment files', required=True)
        parser.add_argument('-f', '--fasta', nargs='+', type=argparse.FileType('r'), help='read inputs for alignment', required=False)
        parser.add_argument('-r', '--read_info', nargs='+', type=str, help='read summary file', required=False)

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

        if not fileExists(args.read_info):
            raise PSToolException("Read info file does not exist: " + str(args.read_info))

        allReadData = DataFrame.parseFromFile(args.read_info, cDelim='\t')

        self.readInfo = {}

        readNameIdx = allReadData.getColumnIndex("READ_NAME")
        readTypeIdx = allReadData.getColumnIndex("TYPE")
        readRunIdx = allReadData.getColumnIndex("RUN_NAME")

        for elem in allReadData.vElements:

            readName = elem[readNameIdx]
            readType = elem[readTypeIdx]
            readRun  = elem[readRunIdx]

            self.readInfo[readName] = (readType, readRun)

    def readSequences(self, args):

        if not fileExists(args.fasta):
            PSToolException("fasta file does not exist: " + args.fasta)

        self.fasta = {}

        for seq in HTSeq.FastaReader( args.fasta ):
            self.fasta[ seq.name ] = seq

    def readGenomeAnnotation(self, args):

        if not fileExists(args.fasta):
            PSToolException("fasta file does not exist: " + args.fasta)

        self.genome_annotation = {}

        for seq in HTSeq.FastaReader( args.fasta ):
            self.fasta[ seq.name ] = seq

    def _makePropDict(self):

        props = {}

        for cigarOp in self.cigar_ops:
            props[cigarOp] = []

        return props

    def prepareInputs(self, args):

        self.readReadInfo(args)
        self.readSequences(args)


        return [x for x in args.sam]

    def execParallel(self, data, environment):



        for readAlignment in HTSeq.SAM_Reader( environment.sam ):

            if readAlignment.aligned:

                readName = readAlignment.read.name
                readSeq = readAlignment.read.seq
                readInterval = readAlignment.iv
                cigarOfRead = readAlignment.cigar

                cigarOp2Length = self._makePropDict()

                for cigarOp in cigarOfRead:

                    if cigarOp.type == 'M':

                        fastaReq = self.fasta[ readInterval.chrom ]
                        fastaSeq = fastaReq[ readInterval.start, readInterval.end ]

                        if readInterval.strand == '-':
                            fastaSeq = fastaSeq.get_reverse_complement()

                        editDistance = self.calculateEditDistance( readSeq, fastaSeq )

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



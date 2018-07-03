import argparse

import pickle
from collections import defaultdict

from porestat.tools.kmer_coverage import KmerHistogram
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

    def __init__(self, parser, subparsers, which):

        super(AlignmentStatisticAnalysisFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):

        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-s', '--sam', nargs='+', type=argparse.FileType('r'), help='alignment files', required=True)
        parser.add_argument('-f', '--fasta', dest='fasta_files', nargs='+', type=argparse.FileType('r'), help='read inputs for alignment', required=True)
        parser.add_argument('-r', '--read-info', nargs='+', type=argparse.FileType('r'), help='read summary file', required=False)

        parser.add_argument('--mc', dest="mc", nargs=1, type=int, default=10, required=False)
        parser.add_argument('--errork', dest="errork", nargs=1, type=int, default=5, required=False)
        parser.add_argument('--perfectk', dest="perfectk", nargs=1, type=int, default=21, required=False)

        parser.add_argument('-q', '--read-type', nargs='+', dest='read_type', action=Fast5TYPEAction, default=None)
        parser.add_argument('-v', '--violin', dest='violin', action='store_true', default=False)

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return AlignmentStatisticAnalysis(simArgs)


"""

--save-parallel-result
/home/mjoppich/dev/data/minion_basecalled/170329_2d_sequencing_run_1.1_p12_pooled/align_summary.pickle
--load-parallel-result
/home/mjoppich/dev/data/minion_basecalled/170329_2d_sequencing_run_1.1_p12_pooled/align_summary.pickle

"""



class AlignmentStatisticAnalysis(ParallelAlignmentPSTReportableInterface):

    def __init__(self, args):

        super(AlignmentStatisticAnalysis, self).__init__(args)

        self.cigar_ops = ['M', 'X', '=', 'D', 'I', 'S', 'X', 'H', 'N', 'P']

        self.readInfoIdx = {}


    def readReadInfo(self, args):

        if args.read_info == None:
            raise PSToolException("Read infos not given ... You must run poreSTAT info first")

        for readInfoFile in args.read_info:

            if not fileExists(readInfoFile.name):
                raise PSToolException("Read info file does not exist: " + str(readInfoFile))

        args.readInfo = None

        for readInfoFile in args.read_info:
            allReadData = DataFrame.parseFromFile(readInfoFile.name, cDelim='\t')

            samIdx = allReadData.addColumn('SAM_READ_NAME')
            readNameIdx = allReadData.getColumnIndex('READ_NAME')
            allReadData.applyByRow(samIdx, lambda x: x[readNameIdx].split(' ')[0])

            if args.readInfo == None:
                args.readInfo = allReadData
            else:
                args.readInfo.merge( allReadData )

    def readSequences(self, args):

        for fastaFile in args.fasta_files:

            if not fileExists(fastaFile.name):
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
            if not fileExists(x.name):
                PSToolException("sam file does not exist: " + str(x))

        return [x for x in args.sam]

    def _createLocalEnvironment(self):
        return { 'BASE_COMPOSITION': defaultdict(Counter),
                 'UNALIGNED_READS': set(),
                 'READ_STATS': {},
                 'KMER_STATS_ALIGNED': {},
                 'KMER_PERF_ALIGNED': {},
                 'READ_INFO': {}
                 }


    def calc_gc(self, readSeq):

        gcCount = readSeq.count("G")
        gcCount += readSeq.count("C")

        return gcCount / len(readSeq)


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
            alignedKmers = defaultdict(lambda: Counter())
            perfKmers = Counter()

            readIdentityCount = 0
            alignmentIdentityCount = 0
            alignmentLength = readAlignment.iv.end-readAlignment.iv.start

            matchedGCContents = []
            matchedLongestLength = None

            fastaReq = globalEnv.fasta[readInterval.chrom]
            fastaReq = fastaReq.toHTSeq()


            for cigarOp in cigarOfRead:

                readSeq = readAlignment.read_as_aligned[cigarOp.query_from:cigarOp.query_to]

                for base in readSeq.seq:
                    localEnv['BASE_COMPOSITION'][cigarOp.type][chr(base)] += 1

                if cigarOp.type == 'M':

                    fastaSeq = fastaReq[ cigarOp.ref_iv.start:cigarOp.ref_iv.end ]

                    editDistance = self.calculateEditDistance( readSeq, fastaSeq )


                    gc_content = self.calc_gc(str(readSeq.seq))
                    matchedGCContents.append( (gc_content, len(readSeq)) )


                    readIdentityCount += len(readSeq) - editDistance
                    alignmentIdentityCount += len(readSeq) - editDistance

                    if len(readSeq) != len(fastaSeq):
                        print(readSeq, fastaSeq, "strange lengths")
                    else:

                        fseq = str(readSeq)
                        aseq = str(fastaSeq)

                        equals = 0
                        mismatches = 0
                        for i in range(0, len(fseq)):

                            if fseq[i] == aseq[i]:
                                equals += 1

                                if mismatches > 0:
                                    cigarOp2Length['Xd'].append(mismatches)
                                    mismatches = 0

                            else:
                                mismatches += 1

                                if equals > 0:
                                    cigarOp2Length['Md'].append(equals)

                                    if i -5 >= 0:
                                        kmerBeforeMM = fseq[i-5:i]
                                        alignedKmers['Md'][kmerBeforeMM] += 1
                                        
                                        
                                    if equals >= 21:
                                        kmerseq = fseq[i-equals:i]
                                        foundkmers = KmerHistogram.calcKmers(kmerseq, 21)
                                        perfKmers += foundkmers



                                    equals = 0

                        if equals > 0:
                            cigarOp2Length['Md'].append(equals)

                        if mismatches > 0:
                            cigarOp2Length['Xd'].append(mismatches)

                    cigarOp2Length['M'].append(cigarOp.size)
                    cigarOp2Length['X'].append(editDistance)
                    cigarOp2Length['='].append(cigarOp.size - editDistance)

                else:

                    kmerStart = cigarOp.query_from - 6
                    kmerEnd = cigarOp.query_from -1

                    if kmerStart >= 0:
                        kmerBeforeError = str(readAlignment.read_as_aligned[kmerStart:kmerEnd])
                        alignedKmers[cigarOp.type][kmerBeforeError] += 1

                    cigarOp2Length[ cigarOp.type ].append(cigarOp.size)

            localEnv['READ_STATS'][(readID, readType)] = cigarOp2Length
            localEnv['KMER_STATS_ALIGNED'][(readID, readType)] = alignedKmers
            localEnv['KMER_PERF_ALIGNED'][(readID, readType)] = perfKmers

            localEnv['READ_INFO'][(readID, readType)] = {
                'READ_IDENTITY': (readIdentityCount /len(readAlignment.read), len(readAlignment.read)),
                'ALIGNMENT_IDENTITY': (alignmentIdentityCount/ alignmentLength, alignmentLength),
                'MATCHED_GC_CONTENT': matchedGCContents,
                'MATCHED_LONGEST_LENGTH': [matchedLongestLength]
            }

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


        for file, data in newResult:

            if file in existResult:
                existResult[file] = mergeDicts(existResult[file], data)
            else:
                existResult[file] = data

        return existResult


    def makeResults(self, allFilesResult, oEnvironment, args):

        for fileName in allFilesResult:

            parallelResult = allFilesResult[fileName]

            readCIGARs = parallelResult['READ_STATS']
            readInfo = parallelResult['READ_INFO']
            cigarKMERs = parallelResult['KMER_STATS_ALIGNED']
            perfectKMERs = parallelResult['KMER_PERF_ALIGNED']

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

            PorePlot.plotBars(overviewByAligned, fileName + ": Alignment Result By Type", "", "Count", pltcfg=args.pltcfg)

            if not self.hasArgument('read_type', args) or not args.read_type:

                totalCounter = defaultdict(list)
                kmerCounter = defaultdict(lambda: Counter())
                perfectCounter = Counter()

                allReadInfos = defaultdict(list)

                for counterPair in readCIGARs:

                    readCounter = readCIGARs[counterPair]
                    readKmers = cigarKMERs[counterPair]

                    allReadInfo = readInfo[counterPair]

                    """
                    {
                        'READ_IDENTITY': (readIdentityCount / len(readAlignment.read), len(readAlignment.read)),
                        'ALIGNMENT_IDENTITY': (alignmentIdentityCount / alignmentLength, alignmentLength),
                        'MATCHED_GC_CONTENT': matchedGCContents,
                        'MATCHED_LONGEST_LENGTH': [matchedLongestLength]
                    }
                    """

                    allReadInfos['READ_IDENTITY'].append(allReadInfo['READ_IDENTITY'])
                    allReadInfos['ALIGNMENT_IDENTITY'].append(allReadInfo['ALIGNMENT_IDENTITY'])
                    allReadInfos['MATCHED_GC_CONTENT'] += allReadInfo['MATCHED_GC_CONTENT']
                    allReadInfos['MATCHED_LONGEST_LENGTH'] += allReadInfo['MATCHED_LONGEST_LENGTH']

                    perfectCounter += perfectKMERs[counterPair]

                    for x in readCounter:
                        totalCounter[x] += readCounter[x]

                    for x in readKmers:
                        kmerCounter[x] += readKmers[x]

                obsInfo = {
                    'RUNID': fileName + "_perfect",
                    'USER_RUN_NAME': "perfect",
                    'KMERCOUNTS': perfectCounter
                }

                kmerDF = KmerHistogram.dfSummary(obsInfo, args.mc)
                kmerDF.setTitle("k-mer histogram: perfectly aligned (k={})".format(args.perfectk))
                args.pltcfg.makeTable(kmerDF)


                for evtype in kmerCounter:

                    obsInfo = {
                        'RUNID': fileName + "_" + str(evtype),
                        'USER_RUN_NAME': evtype,
                        'KMERCOUNTS': kmerCounter[evtype]
                    }

                    kmerDF = KmerHistogram.dfSummary(obsInfo, args.mc)
                    kmerDF.setTitle("k-mer histogram: sequence before CIGAR " + evtype + " (k={})".format(args.errork))
                    args.pltcfg.makeTable(kmerDF)




                """
                READ IDENTITY, ALIGNMENT IDENTITY, MATCHED GC CONTENT, MAX MATCHED PART
                """

                xReadID = [x[0] for x in allReadInfos['READ_IDENTITY']]
                yReadID = [x[1] for x in allReadInfos['READ_IDENTITY']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Identity vs Read Length", "Read Identity", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos['ALIGNMENT_IDENTITY']]
                yReadID = [x[1] for x in allReadInfos['ALIGNMENT_IDENTITY']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Identity vs Alignment Length", "Alignment Identity", "Alignment Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos['MATCHED_GC_CONTENT']]
                yReadID = [x[1] for x in allReadInfos['MATCHED_GC_CONTENT']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "GC Content vs Aligned Seq Length", "GC Content", "Aligned Sequence Length [bp]", pltcfg=args.pltcfg)

                if self.hasArgument('violin', args) and args.violin:
                    PorePlot.plotViolin(totalCounter, sorted([x for x in totalCounter]), "CIGARs: all types", pltcfg=args.pltcfg, shareX = False, shareY=False)
                else:
                    PorePlot.plotBoxplot(totalCounter, sorted([x for x in totalCounter]), "CIGARs: all types", pltcfg=args.pltcfg)


            return


            totalCounter = defaultdict(lambda: defaultdict(list))
            kmersByReadType = defaultdict(lambda: defaultdict(lambda: Counter()))
            perfectCounter = defaultdict(lambda: Counter())

            allReadInfos = defaultdict(lambda: defaultdict(list))

            for counterPair in readCIGARs:

                readID = counterPair[0]
                readType = counterPair[1]
                readCounter = readCIGARs[counterPair]
                readKmers = cigarKMERs[counterPair]

                allReadInfo = readInfo[counterPair]
                allReadInfos[readType] = mergeDicts(allReadInfos[readType], allReadInfo)

                allReadInfos[readType]['READ_IDENTITY'].append(allReadInfo['READ_IDENTITY'])
                allReadInfos[readType]['ALIGNMENT_IDENTITY'].append(allReadInfo['ALIGNMENT_IDENTITY'])
                allReadInfos[readType]['MATCHED_GC_CONTENT'] += allReadInfo['MATCHED_GC_CONTENT']
                allReadInfos[readType]['MATCHED_LONGEST_LENGTH'] += allReadInfo['MATCHED_LONGEST_LENGTH']

                perfectCounter[readType] += perfectKMERs[counterPair]

                for x in readCounter:
                    totalCounter[readType][x] += readCounter[x]

                for x in readKmers:
                    kmersByReadType[readType][x] += readKmers[x]


            for readType in totalCounter:

                print(readType)
                print(self.hasArgument('violin', args),args.violin)

                obsInfo = {
                    'RUNID': fileName + "_" + str(readType) + "_perfect",
                    'USER_RUN_NAME': "perfect",
                    'KMERCOUNTS': perfectCounter[readType]
                }
                kmerDF = KmerHistogram.dfSummary(obsInfo, args.mc)
                kmerDF.setTitle("k-mer histogram: perfectly aligned (k={}, {})".format(args.perfectk, str(readType)))
                args.pltcfg.makeTable(kmerDF)

                for rtype in kmersByReadType[readType]:

                    print(rtype)

                    obsInfo = {
                        'RUNID': fileName + "_" + str(readType) + "_" + str(rtype),
                        'USER_RUN_NAME': rtype,
                        'KMERCOUNTS': kmersByReadType[readType][rtype]
                    }
                    kmerDF = KmerHistogram.dfSummary(obsInfo, args.mc)
                    kmerDF.setTitle("k-mer histogram: sequence before CIGAR " + rtype + " (k={}, {})".format(args.errork, str(readType)))
                    args.pltcfg.makeTable(kmerDF)






                """
                READ IDENTITY, ALIGNMENT IDENTITY, MATCHED GC CONTENT, MAX MATCHED PART
                """

                xReadID = [x[0] for x in allReadInfos[readType]['READ_IDENTITY']]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_IDENTITY']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Identity vs Read Length ({})".format(str(readType)), "Read Identity", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['ALIGNMENT_IDENTITY']]
                yReadID = [x[1] for x in allReadInfos[readType]['ALIGNMENT_IDENTITY']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Identity vs Alignment Length  ({})".format(str(readType)), "Alignment Identity", "Alignment Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['MATCHED_GC_CONTENT']]
                yReadID = [x[1] for x in allReadInfos[readType]['MATCHED_GC_CONTENT']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "GC Content vs Aligned Seq Length  ({})".format(str(readType)), "GC Content", "Aligned Sequence Length [bp]", pltcfg=args.pltcfg)


                if self.hasArgument('violin', args) and args.violin:
                    PorePlot.plotViolin(totalCounter[readType], sorted([x for x in totalCounter[readType]]), "CIGARs: " + str(readType),
                                        pltcfg=args.pltcfg, shareX = False, shareY=False)
                else:
                    PorePlot.plotBoxplot(totalCounter[readType], sorted([x for x in totalCounter[readType]]), "CIGARs: " + str(readType),
                                         pltcfg=args.pltcfg)
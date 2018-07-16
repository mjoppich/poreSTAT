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
from ..utils.DataFrame import DataFrame, DataRow

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
                 'READ_INFO': {},
                 "NT_SUBST": Counter()
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

        alignSubst = Counter()

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
            refMatchedGCContents = []
            matchedLongestLength = None

            fastaReq = globalEnv.fasta[readInterval.chrom]
            fastaReq = fastaReq.toHTSeq()


            readLength = len(readAlignment.read)
            alignQual = readAlignment.aQual
            seqQual = (min(readAlignment.read.qual), self.mean(readAlignment.read.qual), max(readAlignment.read.qual))
            refLength = alignmentLength

            for cigarOp in cigarOfRead:

                readSeq = readAlignment.read_as_aligned[cigarOp.query_from:cigarOp.query_to]

                for base in readSeq.seq:
                    localEnv['BASE_COMPOSITION'][cigarOp.type][chr(base)] += 1

                if cigarOp.type == 'M':

                    fastaSeq = fastaReq[ cigarOp.ref_iv.start:cigarOp.ref_iv.end ]

                    editDistance = self.calculateEditDistance( readSeq, fastaSeq )


                    gc_content = self.calc_gc(str(readSeq.seq))
                    matchedGCContents.append( (gc_content, len(readSeq)) )

                    ref_gc_content = self.calc_gc(str(fastaSeq))
                    refMatchedGCContents.append((ref_gc_content, len(fastaSeq)))


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

                                #mismatch
                                refNT = aseq[i]
                                readNT = fseq[i]

                                alignSubst[(readNT, refNT)] += 1

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

            localEnv['NT_SUBST'] += alignSubst

            localEnv['READ_INFO'][(readID, readType)] = {
                'READ_IDENTITY': (readIdentityCount /len(readAlignment.read), len(readAlignment.read)),
                'ALIGNMENT_IDENTITY': (alignmentIdentityCount/ alignmentLength, alignmentLength),
                'READ_MATCHED_GC_CONTENT': matchedGCContents,
                'REF_MATCHED_GC_CONTENT': refMatchedGCContents,
                'MATCHED_LONGEST_LENGTH': [matchedLongestLength],

                'READ_LENGTH': readLength,
                'REF_LENGTH': refLength,
                'SEQ_QUAL': seqQual,
                'ALIGN_QUAL': alignQual,
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


    @classmethod
    def mean(cls, elems):

        return sum(elems) / len(elems)

    def makeResults(self, allFilesResult, oEnvironment, args):

        for fileName in allFilesResult:

            parallelResult = allFilesResult[fileName]

            readCIGARs = parallelResult['READ_STATS']
            readInfo = parallelResult['READ_INFO']
            cigarKMERs = parallelResult['KMER_STATS_ALIGNED']
            perfectKMERs = parallelResult['KMER_PERF_ALIGNED']

            ntSubstitutions = parallelResult['NT_SUBST']

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


            totalCounter = defaultdict(list)
            kmerCounter = defaultdict(lambda: Counter())
            perfectCounter = Counter()

            allReadInfos = defaultdict(list)
            totalReadBases = 0

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
                allReadInfos['MATCHED_LONGEST_LENGTH'] += allReadInfo['MATCHED_LONGEST_LENGTH']

                allReadInfos['READ_MATCHED_GC_CONTENT'] += allReadInfo['READ_MATCHED_GC_CONTENT']
                allReadInfos['REF_MATCHED_GC_CONTENT'] += allReadInfo['REF_MATCHED_GC_CONTENT']


                allReadInfos['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH'].append( ( sum([x[0]*x[1] for x in allReadInfo['READ_MATCHED_GC_CONTENT']]) / sum([x[1] for x in allReadInfo['READ_MATCHED_GC_CONTENT']])  , allReadInfo['READ_LENGTH'])  )
                allReadInfos['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH'].append(  (  sum([x[0] * x[1] for x in allReadInfo['REF_MATCHED_GC_CONTENT']]) / sum([x[1] for x in allReadInfo['REF_MATCHED_GC_CONTENT']]), allReadInfo['REF_LENGTH'])    )

                totalReadBases += allReadInfo['READ_LENGTH']

                allReadInfos['READ_ALIGN_QUAL_TO_READ_LENGTH'].append(  (allReadInfo['ALIGN_QUAL'],  allReadInfo['READ_LENGTH']) )
                allReadInfos['SEQ_QUAL_MIN_TO_READ_LENGTH'].append(     (allReadInfo['SEQ_QUAL'][0], allReadInfo['READ_LENGTH']) )
                allReadInfos['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN'].append( (allReadInfo['ALIGN_QUAL'],  allReadInfo['SEQ_QUAL'][0]) )
                allReadInfos['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN'].append(
                    (allReadInfo['ALIGN_QUAL'], allReadInfo['SEQ_QUAL'][1]))

                perfectCounter += perfectKMERs[counterPair]

                for x in readCounter:
                    totalCounter[x] += readCounter[x]

                for x in readKmers:
                    kmerCounter[x] += readKmers[x]






            overallIdentity = (sum(totalCounter["="]) / totalReadBases) * 100
            alignedIdentity = (sum(totalCounter["="]) / (sum(totalCounter['=']) + sum(totalCounter['X']))) * 100

            alignedIdentity = (sum(totalCounter["="]) / sum(totalCounter['M'])) * 100
            insertedPerAligned = (sum(totalCounter["I"]) / sum(totalCounter['M'])) * 100
            deletedPerAligned = (sum(totalCounter["D"]) / sum(totalCounter['M'])) * 100
            substitutedPerAligned = (sum(totalCounter["X"]) / sum(totalCounter['M'])) * 100


            fmtStr = "{:.4f}"

            statdf = DataFrame()
            statdf.setTitle("Alignment Statistics")
            statdf.addColumns(['Name', 'Value'])

            statdf.addRow(DataRow.fromDict({
                'Name': 'Overall Identity',
                'Value': fmtStr.format(overallIdentity)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Aligned Identity',
                'Value': fmtStr.format(alignedIdentity)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Identical bases per 100 aligned bases', # = / M * 100
                'Value': fmtStr.format(alignedIdentity)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Inserted bases per 100 aligned bases',
                'Value': fmtStr.format(insertedPerAligned)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Deleted bases per 100 aligned bases',
                'Value': fmtStr.format(deletedPerAligned)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Substituted bases per 100 aligned bases',
                'Value': fmtStr.format(substitutedPerAligned)
            }))

            addedGC = sum([x[0]*x[1] for x in allReadInfos['READ_MATCHED_GC_CONTENT']])
            addedLength = sum([x[1] for x in allReadInfos['READ_MATCHED_GC_CONTENT']])

            statdf.addRow(DataRow.fromDict({
                'Name': 'Avg GC content of Aligned Read Sequence',
                'Value': fmtStr.format(addedGC/addedLength)
            }))

            gcContentRef = 0
            refLength = 0

            for chrom in oEnvironment.fasta:

                fastaReq = oEnvironment.fasta[chrom]
                fastaReq = fastaReq.toHTSeq()

                strSeq = str(fastaReq.seq)
                gcContentRef += self.calc_gc(strSeq) * len(strSeq)
                refLength += len(strSeq)

            statdf.addRow(DataRow.fromDict({
                'Name': 'Avg GC content of Reference Sequence',
                'Value': fmtStr.format(gcContentRef/refLength)
            }))


            args.pltcfg.makeTable(statdf)



            allSubst = sum([ntSubstitutions[x] for x in ntSubstitutions])

            substDF = DataFrame()
            substDF.setTitle("Substitution Statistics (read -> ref)")
            substDF.addColumns(['Substitution', 'Abs Count', 'Rel Subst'])

            for subst in ntSubstitutions:
                substDF.addRow(DataRow.fromDict({
                    'Substitution': subst[0]+' -> ' + subst[1],
                    'Abs Count': fmtStr.format(ntSubstitutions[subst]),
                    'Rel Subst': fmtStr.format(ntSubstitutions[subst]/allSubst)
                }))

            args.pltcfg.makeTable(substDF)



            if not self.hasArgument('read_type', args) or not args.read_type:

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

                xReadID = [x[0] for x in allReadInfos['READ_MATCHED_GC_CONTENT']]
                yReadID = [x[1] for x in allReadInfos['READ_MATCHED_GC_CONTENT']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Matched GC Content vs Aligned Seq Length", "GC Content", "Aligned Sequence Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos['REF_MATCHED_GC_CONTENT']]
                yReadID = [x[1] for x in allReadInfos['REF_MATCHED_GC_CONTENT']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Reference Matched GC Content vs Aligned Seq Length", "GC Content", "Aligned Sequence Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH']]
                yReadID = [x[1] for x in allReadInfos['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Matched GC Content vs Read Length", "GC Content", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH']]
                yReadID = [x[1] for x in allReadInfos['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Reference Matched GC Content vs Read Length", "GC Content", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos['READ_ALIGN_QUAL_TO_READ_LENGTH']]
                yReadID = [x[1] for x in allReadInfos['READ_ALIGN_QUAL_TO_READ_LENGTH']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Align Quality vs Read Length", "Read Align Quality", "Read Lenght [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos['SEQ_QUAL_MIN_TO_READ_LENGTH']]
                yReadID = [x[1] for x in allReadInfos['SEQ_QUAL_MIN_TO_READ_LENGTH']]

                print("min read length", min(yReadID), "max", max(yReadID))
                print("xvals", xReadID[yReadID.index(min(yReadID))], xReadID[yReadID.index(max(yReadID))])

                PorePlot.plot_scatter_densities(xReadID, yReadID, "Minimum Sequence Quality vs Read Length (min value of read)", "Minimum Sequence Quality", "Read Length [bp]", pltcfg=args.pltcfg)


                xReadID = [x[0] for x in allReadInfos['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN']]
                yReadID = [x[1] for x in allReadInfos['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Quality vs Minimum Sequence Quality", "Alignment Quality", "Sequence Quality (min value of read)", pltcfg=args.pltcfg)


                xReadID = [x[0] for x in allReadInfos['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN']]
                yReadID = [x[1] for x in allReadInfos['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Quality vs Mean Sequence Quality", "Alignment Quality", "Sequence Quality (mean value of read)", pltcfg=args.pltcfg)


                if self.hasArgument('violin', args) and args.violin:
                    PorePlot.plotViolin(totalCounter, sorted([x for x in totalCounter]), "CIGARs: all types", pltcfg=args.pltcfg, shareX = False, shareY=False)
                else:
                    PorePlot.plotBoxplot(totalCounter, sorted([x for x in totalCounter]), "CIGARs: all types", pltcfg=args.pltcfg)


            totalCounter = defaultdict(lambda: defaultdict(list))
            kmersByReadType = defaultdict(lambda: defaultdict(lambda: Counter()))
            perfectCounter = defaultdict(lambda: Counter())

            allReadInfos = defaultdict(lambda: defaultdict(list))

            totalReadBases = Counter()

            for counterPair in readCIGARs:

                readID = counterPair[0]
                readType = counterPair[1]
                readCounter = readCIGARs[counterPair]
                readKmers = cigarKMERs[counterPair]

                allReadInfo = readInfo[counterPair]
                #allReadInfos[readType] = mergeDicts(allReadInfos[readType], allReadInfo)

                allReadInfos[readType]['READ_IDENTITY'].append(allReadInfo['READ_IDENTITY'])
                allReadInfos[readType]['ALIGNMENT_IDENTITY'].append(allReadInfo['ALIGNMENT_IDENTITY'])

                totalReadBases[readType] += allReadInfo['READ_LENGTH']

                allReadInfos[readType]['READ_MATCHED_GC_CONTENT'] += allReadInfo['READ_MATCHED_GC_CONTENT']
                allReadInfos[readType]['REF_MATCHED_GC_CONTENT'] += allReadInfo['REF_MATCHED_GC_CONTENT']

                allReadInfos[readType]['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH'].append((sum(
                    [x[0] * x[1] for x in allReadInfo['READ_MATCHED_GC_CONTENT']]) / sum(
                    [x[1] for x in allReadInfo['READ_MATCHED_GC_CONTENT']]), allReadInfo['READ_LENGTH']))
                allReadInfos[readType]['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH'].append((sum(
                    [x[0] * x[1] for x in allReadInfo['REF_MATCHED_GC_CONTENT']]) / sum(
                    [x[1] for x in allReadInfo['REF_MATCHED_GC_CONTENT']]), allReadInfo['REF_LENGTH']))

                allReadInfos[readType]['READ_ALIGN_QUAL_TO_READ_LENGTH'].append(
                    (allReadInfo['ALIGN_QUAL'], allReadInfo['READ_LENGTH']))
                allReadInfos[readType]['SEQ_QUAL_MIN_TO_READ_LENGTH'].append(
                    (allReadInfo['SEQ_QUAL'][0], allReadInfo['READ_LENGTH']))
                allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN'].append(
                    (allReadInfo['ALIGN_QUAL'], allReadInfo['SEQ_QUAL'][0]))

                allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN'].append(
                    (allReadInfo['ALIGN_QUAL'], allReadInfo['SEQ_QUAL'][1]))
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

                overallIdentity = (sum(totalCounter[readType]["="]) / totalReadBases[readType]) * 100
                alignedIdentity = (sum(totalCounter[readType]["="]) / (sum(totalCounter[readType]['=']) + sum(totalCounter[readType]['X']))) * 100

                alignedIdentity = (sum(totalCounter[readType]["="]) / sum(totalCounter[readType]['M'])) * 100
                insertedPerAligned = (sum(totalCounter[readType]["I"]) / sum(totalCounter[readType]['M'])) * 100
                deletedPerAligned = (sum(totalCounter[readType]["D"]) / sum(totalCounter[readType]['M'])) * 100
                substitutedPerAligned = (sum(totalCounter[readType]["X"]) / sum(totalCounter[readType]['M'])) * 100

                fmtStr = "{:.4f}"

                statdf = DataFrame()
                statdf.setTitle("Alignment Statistics ({})".format(readType))
                statdf.addColumns(['Name', 'Value'])

                statdf.addRow(DataRow.fromDict({
                    'Name': 'Overall Identity',
                    'Value': fmtStr.format(overallIdentity)
                }))

                statdf.addRow(DataRow.fromDict({
                    'Name': 'Aligned Identity',
                    'Value': fmtStr.format(alignedIdentity)
                }))

                statdf.addRow(DataRow.fromDict({
                    'Name': 'Identical bases per 100 aligned bases',  # = / M * 100
                    'Value': fmtStr.format(alignedIdentity)
                }))

                statdf.addRow(DataRow.fromDict({
                    'Name': 'Inserted bases per 100 aligned bases',
                    'Value': fmtStr.format(insertedPerAligned)
                }))

                statdf.addRow(DataRow.fromDict({
                    'Name': 'Deleted bases per 100 aligned bases',
                    'Value': fmtStr.format(deletedPerAligned)
                }))

                statdf.addRow(DataRow.fromDict({
                    'Name': 'Substitutions bases per 100 aligned bases',
                    'Value': fmtStr.format(substitutedPerAligned)
                }))

                args.pltcfg.makeTable(statdf)





                """
                READ IDENTITY, ALIGNMENT IDENTITY, MATCHED GC CONTENT, MAX MATCHED PART
                """

                xReadID = [x[0] for x in allReadInfos[readType]['READ_IDENTITY']]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_IDENTITY']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Identity vs Read Length ({})".format(str(readType)), "Read Identity", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['ALIGNMENT_IDENTITY']]
                yReadID = [x[1] for x in allReadInfos[readType]['ALIGNMENT_IDENTITY']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Identity vs Alignment Length  ({})".format(str(readType)), "Alignment Identity", "Alignment Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['READ_MATCHED_GC_CONTENT']]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_MATCHED_GC_CONTENT']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Matched GC Content vs Aligned Seq Length ({})".format(str(readType)), "GC Content", "Aligned Sequence Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['REF_MATCHED_GC_CONTENT']]
                yReadID = [x[1] for x in allReadInfos[readType]['REF_MATCHED_GC_CONTENT']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Reference Matched GC Content vs Aligned Seq Length ({})".format(str(readType)), "GC Content", "Aligned Sequence Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH']]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Matched GC Content vs Read Length ({})".format(str(readType)), "GC Content", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH']]
                yReadID = [x[1] for x in allReadInfos[readType]['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Reference Matched GC Content vs Read Length ({})".format(str(readType)), "GC Content", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_READ_LENGTH']]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_READ_LENGTH']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Align Quality vs Read Length ({})".format(str(readType)), "Read Align Quality", "Read Lenght [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['SEQ_QUAL_MIN_TO_READ_LENGTH']]
                yReadID = [x[1] for x in allReadInfos[readType]['SEQ_QUAL_MIN_TO_READ_LENGTH']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Minimum Sequence Quality vs Read Length (min value of read, {})".format(str(readType)), "Minimum Sequence Quality", "Read Length [bp]", pltcfg=args.pltcfg)


                xReadID = [x[0] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN']]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Quality vs Minimum Sequence Quality ({})".format(str(readType)), "Alignment Quality", "Sequence Quality (min value of read)", pltcfg=args.pltcfg)


                xReadID = [x[0] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN']]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Quality vs Mean Sequence Quality", "Alignment Quality", "Sequence Quality (mean value of read)", pltcfg=args.pltcfg)



                if self.hasArgument('violin', args) and args.violin:
                    PorePlot.plotViolin(totalCounter[readType], sorted([x for x in totalCounter[readType]]), "CIGARs: " + str(readType),
                                        pltcfg=args.pltcfg, shareX = False, shareY=False)
                else:
                    PorePlot.plotBoxplot(totalCounter[readType], sorted([x for x in totalCounter[readType]]), "CIGARs: " + str(readType),
                                         pltcfg=args.pltcfg)
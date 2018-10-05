import HTSeq
import argparse

import pickle
import sys
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

import ctypes
from collections import Counter

lib = ctypes.cdll.LoadLibrary('../clib/lib/libfoo.so')


class COUNTER_PAIR(ctypes.Structure):
    _fields_ = [
        ("key", ctypes.c_char_p),
        ("value", ctypes.c_uint32)
    ]


class SEQ_COVERAGE(ctypes.Structure):
    _fields_ = [
        ("seq", ctypes.c_char_p),
        ("start", ctypes.c_uint32),
        ("end", ctypes.c_uint32)
    ]


class NT_SUBSTITUTION(ctypes.Structure):
    _fields_ = [
        ("src", ctypes.c_char),
        ("tgt", ctypes.c_char),
        ("count", ctypes.c_uint32)
    ]


class READ_ALIGNMENT(ctypes.Structure):
    _fields_ = [
        ("pReadID", ctypes.c_char_p),
        ("aligned", ctypes.c_bool),
        ("iReadLength", ctypes.c_uint32),
        ("iRefLength", ctypes.c_uint32),
        ("iAlignQual", ctypes.c_uint32),
        ("fSeqQual_min", ctypes.c_float),
        ("fSeqQual_median", ctypes.c_float),
        ("fSeqQual_max", ctypes.c_float),
        ("iLongestMatched", ctypes.c_uint32),
        ("fReadGCContent", ctypes.c_float),
        ("fRefGCContent", ctypes.c_float),
        ("fReadIdentity", ctypes.c_float),
        ("fRefIdentity", ctypes.c_float),

        ("PERF_KMER_COUNT", ctypes.c_uint32),
        ("PERF_KMERS", ctypes.c_void_p),

        ("CIGAR_COUNT", ctypes.c_uint32),
        ("CIGARS", ctypes.c_char_p),
        ("CIGAR_VEC_LENGTHS", ctypes.POINTER(ctypes.c_uint32)),
        ("CIGAR_LENGTHS", ctypes.POINTER(ctypes.POINTER(ctypes.c_uint32))),

        ("MM2KMER_VEC_LENGTHS", ctypes.POINTER(ctypes.c_uint32)),
        ("MM2KMER_LENGTHS", ctypes.POINTER(ctypes.POINTER(ctypes.c_char_p))),

        ("COV_COUNT", ctypes.c_uint32),
        ("COVERAGES", ctypes.c_void_p),

        ("NT_SUBST_COUNT", ctypes.c_uint32),
        ("NT_SUBST", ctypes.c_void_p),

    ]


class ReadAlignmentStats:

    def __init__(self):

        self.read_id = None
        self.aligned = None
        self.read_type = None

        self.read_length = None
        self.ref_length = None
        self.align_quality = None
        self.seqqual_min = None
        self.seqqual_median = None
        self.seqqual_max = None
        self.longest_matched_stretch = None
        self.read_gc = None
        self.ref_gc = None
        self.read_identity = None
        self.ref_identity = None
        self.nt_substitutions = None
        self.cigar_lengths = None
        self.perfect_kmers = None
        self.coverages = None
        self.mm2kmer = None

class CAlignmentStatistics(object):
    def __init__(self):
        lib.AlignmentStatistics_new.argtypes = []
        lib.AlignmentStatistics_new.restype = ctypes.c_void_p

        lib.AlignmentStatistics_process.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
        lib.AlignmentStatistics_process.restype = ctypes.c_void_p

        lib.AlignmentStatistics_load_fasta.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        lib.AlignmentStatistics_load_fasta.restype = ctypes.c_void_p

        lib.AlignmentStatistics_ReadStats.argtypes = [ctypes.c_void_p]
        lib.AlignmentStatistics_ReadStats.restype = ctypes.c_void_p

        lib.AlignmentStatistics_ReadStats_size.argtypes = [ctypes.c_void_p]
        lib.AlignmentStatistics_ReadStats_size.restype = ctypes.c_uint32

        self.obj = lib.AlignmentStatistics_new()


    def loadFasta(self, fastaPath):

        fp = fastaPath.encode('utf-8')

        charStr = ctypes.cast(fp, ctypes.c_char_p)

        lib.AlignmentStatistics_load_fasta(self.obj, charStr)

    def processFiles(self, vals):

        vals = [x.encode('utf-8') for x in vals]
        c_array = (ctypes.c_char_p * len(vals))(*vals)

        charPArray = ctypes.cast(c_array, ctypes.POINTER(ctypes.c_char_p))

        print(type(charPArray))

        retSize = lib.AlignmentStatistics_process(self.obj, ctypes.c_int(len(vals)), charPArray)
        print(retSize)

        retVecSize = lib.AlignmentStatistics_ReadStats_size(self.obj)
        print("Vector Size")
        print(retVecSize)

        retVec = lib.AlignmentStatistics_ReadStats(self.obj)

        s = ctypes.cast(retVec, ctypes.POINTER(READ_ALIGNMENT * retVecSize))

        allElems = [x for x in s.contents]

        readResults = {}

        print("Printing vector elems")
        for i in range(0, retVecSize):

            elem = allElems[i]

            """

            CIGAR2LEN
            """
            cigars = [x for x in elem.CIGARS.decode()[0:elem.CIGAR_COUNT]]
            vecLengths = [x for x in ctypes.cast(elem.CIGAR_VEC_LENGTHS,
                                                 ctypes.POINTER(ctypes.c_uint32 * elem.CIGAR_COUNT)).contents]
            cigarVecs = [x for x in
                         ctypes.cast(elem.CIGAR_LENGTHS, ctypes.POINTER(ctypes.c_void_p * elem.CIGAR_COUNT)).contents]

            cigar2len = {}
            for c in range(0, len(cigars)):
                cigarLen = [x for x in
                            ctypes.cast(cigarVecs[c], ctypes.POINTER(ctypes.c_uint32 * vecLengths[c])).contents]

                cigar2len[cigars[c]] = cigarLen

            print("CIGAR2LEN")
            for x in cigar2len:
                print(x, cigar2len[x])

            """

            PERF KMERS
            """
            perfKmers = [x for x in
                         ctypes.cast(elem.PERF_KMERS, ctypes.POINTER(COUNTER_PAIR * elem.PERF_KMER_COUNT)).contents]

            perfKmerCounter = Counter()
            for x in perfKmers:
                perfKmerCounter[x.key.decode()] = x.value

            print("PERF KMER")
            for x in perfKmerCounter:
                print(x, perfKmerCounter[x])

            """

            CIGAR2LEN
            """
            mmLengths = [x for x in ctypes.cast(elem.MM2KMER_VEC_LENGTHS,
                                                ctypes.POINTER(ctypes.c_uint32 * elem.CIGAR_COUNT)).contents]
            mm2kmers = [x for x in
                        ctypes.cast(elem.MM2KMER_LENGTHS, ctypes.POINTER(ctypes.c_void_p * elem.CIGAR_COUNT)).contents]

            mm2kmer = {}
            for c in range(0, len(cigars)):
                mmKmers = [x.decode() for x in
                           ctypes.cast(mm2kmers[c], ctypes.POINTER(ctypes.c_char_p * mmLengths[c])).contents]

                mm2kmer[cigars[c]] = mmKmers

            print("MM2KMER")
            for x in mm2kmer:
                print(x, mm2kmer[x])

            """

            COVERAGES
            """
            covIntervals = [x for x in
                            ctypes.cast(elem.COVERAGES, ctypes.POINTER(SEQ_COVERAGE * elem.COV_COUNT)).contents]

            # TODO add this to global array

            """

            NT SUBSTITUTION
            """
            ntSubst = [x for x in
                       ctypes.cast(elem.NT_SUBST, ctypes.POINTER(NT_SUBSTITUTION * elem.NT_SUBST_COUNT)).contents]

            substCounter = Counter()
            for x in ntSubst:
                substCounter[(x.src.decode(), x.tgt.decode())] = x.count

            print("NT SUBST")
            for x in substCounter:
                print(x, substCounter[x])

            print(allElems[i].pReadID)

            readStat = ReadAlignmentStats()

            readStat.read_id = elem.pReadID.decode()
            readStat.aligned = elem.aligned
            readStat.read_length = elem.iReadLength
            readStat.ref_length = elem.iRefLength
            readStat.align_quality = elem.iAlignQual
            readStat.seqqual_min = elem.fSeqQual_min
            readStat.seqqual_median = elem.fSeqQual_median
            readStat.seqqual_max = elem.fSeqQual_max
            readStat.longest_matched_stretch = elem.iLongestMatched
            readStat.read_gc = elem.fReadGCContent
            readStat.ref_gc = elem.fRefGCContent
            readStat.read_identity = elem.fReadIdentity
            readStat.ref_identity = elem.fRefIdentity
            readStat.nt_substitutions = substCounter
            readStat.cigar_lengths = cigar2len
            readStat.perfect_kmers = perfKmerCounter
            readStat.coverages = covIntervals
            readStat.mm2kmer = mm2kmer

            readResults[readStat.read_id] = readStat


        return readResults



class AlignmentStatisticAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(AlignmentStatisticAnalysisFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):

        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-s', '--sam', nargs='+', type=argparse.FileType('r'), help='alignment files', required=True)
        parser.add_argument('-f', '--fasta', dest='fasta_files', nargs='+', type=argparse.FileType('r'), help='read inputs for alignment', required=True)
        parser.add_argument('-r', '--read-info', nargs='+', type=argparse.FileType('r'), help='read summary file', required=False)

        parser.add_argument('-g', '--gtf', type=argparse.FileType('r'), help="gtf file for features", required=False, default=None)

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

            sys.stderr.write("Read infos not given ... You must run poreSTAT info first")
            #raise PSToolException("Read infos not given ... You must run poreSTAT info first")

        else:

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
        return None


    def calc_gc(self, readSeq):

        gcCount = readSeq.count("G")
        gcCount += readSeq.count("C")

        return gcCount / len(readSeq)


    def getReadInfo(self, readName):

        readInfo = self.args.readInfo.findRow('SAM_READ_NAME', readName)
        if readInfo != None:
            readType = Fast5TYPE[ readInfo[ "TYPE" ] ]
            readID = readInfo[ "READ_ID" ]

        return readType,readID


    def addCoverage(self, chrom, start, end):

        if not chrom in self.genomicArray.chrom_vectors:
            self.genomicArray.add_chrom(chrom)

        self.genomicArray[HTSeq.GenomicInterval( chrom, start, end, "." )] += 1


    def execParallel(self, data, environment):
        """

        :param data: input files
        :param environment: usually NONE
        :return: None
        """
        statMaker = CAlignmentStatistics()

        for fastaFile in self.args.fasta:
            statMaker.loadFasta(fastaFile.name)


        alignResults = {}

        for samFile in data:
            readStats = statMaker.processFiles([samFile])

            for readID in readStats:
                readStats[readID].read_type, _ = self.getReadInfo(readID)

            alignResults[samFile] = readStats


        return alignResults




    @classmethod
    def mean(cls, elems):

        return sum(elems) / len(elems)

    def makeResults(self, allFilesResult, oEnvironment, args):
        """



        :param allFilesResult:
        :param oEnvironment:
        :param args:
        :return:
        """

        print("Preparing results")



        """
        
        FASTA REFERENCE INFOS
        
        """

        gcContentRef = 0
        refLength = 0

        fastaSeqs = {}
        for fastaFile in self.args.fasta:

            for seq in HTSeq.FastaReader(fastaFile):
                fastaSeqs[seq.name] = str(seq.seq)

                strSeq = str(seq.seq)
                gcContentRef += self.calc_gc(strSeq) * len(strSeq)
                refLength += len(strSeq)



        """
        
        GTF/FEATURE INFOS
        
        """

        if args.gtf:

            allFeatureArrays = defaultdict(lambda: HTSeq.GenomicArray( [], stranded=False ) )

            gtf_file = HTSeq.GFF_Reader( args.gtf.name, end_included = True )


            for feature in gtf_file:

                if not feature.type in allFeatureArrays:

                    print("Adding feature", feature.type)

                ga = allFeatureArrays[ feature.type ]

                if not feature.iv.chrom in ga.chrom_vectors:
                    ga.add_chrom(feature.iv.chrom)

                ga[feature.iv] = 1

                allFeatureArrays[feature.type] = ga

            for type in allFeatureArrays:
                print("GTF Feature Type:", type)

        """
        
        CREATING STATS PER samFile
        
        """


        for fileName in allFilesResult:


            args.pltcfg.addHTMLPlot("<h2>" + str(fileName) + "</h2>\n")

            allReads = allFilesResult[fileName]

            alignedIDs = set()
            unalignedIDs = set()

            allAlignedReads = {}
            for x in allReads:

                elem = allReads[x]

                if elem.aligned == True:
                    allAlignedReads[x] = elem
                    alignedIDs.add(elem.read_id)
                else:
                    unalignedIDs.add(elem.read_id)




            parallelResult = allFilesResult[fileName]

            #readCIGARs = parallelResult['READ_STATS']
            #readInfo = parallelResult['READ_INFO']
            #cigarKMERs = parallelResult['KMER_STATS_ALIGNED']
            #perfectKMERs = parallelResult['KMER_PERF_ALIGNED']
            #ntSubstitutions = parallelResult['NT_SUBST']
            #coverageArray = parallelResult['READS_COVERAGE']
            additionalUnalignedReads = set()

            if self.hasArgument('readInfo', self.args) and self.args.readInfo != None:

                for row in self.args.readInfo:

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

            if self.hasArgument('readInfo', self.args) and args.readInfo != None:

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

            if len(overviewByAligned) > 0:
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




            args.pltcfg.addHTMLPlot("<h3>Feature-based Coverage</h3>\n")

            if args.gtf and len(allFeatureArrays) > 0:


                fmtStr = "{:.4f}"

                featureDF = DataFrame()
                featureDF.setTitle("<h3>Alignment Feature Statistics</h3>")
                featureDF.addColumns(['Name', 'Value'])



                for type in allFeatureArrays:

                    featureCovered = 0
                    readsCovered = 0


                    fga = allFeatureArrays[type]

                    for featurePart in fga.steps():

                        partCovered = featurePart[0].length
                        partCount = featurePart[1]

                        if partCount > 0:

                            featureCovered += partCovered

                            #for this part fetch steps in read coverage

                            for readPart in coverageArray[featurePart[0]].steps():

                                readCovered = readPart[0].length
                                readCount = readPart[1]

                                if readCount > 0:
                                    readsCovered += readCovered

                    featureDF.addRow(DataRow.fromDict({
                        'Name': 'Featuretype Coverage ({})'.format(type),
                        'Value': fmtStr.format(featureCovered/refLength)
                    }))

                    featureDF.addRow(DataRow.fromDict({
                        'Name': 'Featuretype Covered By Reads ({})'.format(type),
                        'Value': fmtStr.format(readsCovered/featureCovered)
                    }))

                args.pltcfg.makeTable(featureDF)

            else:

                args.pltcfg.addHTMLPlot("<p>No feature annotation supplied.</p>\n")


            args.pltcfg.addHTMLPlot("<h3>Alignment Statistics</h3>\n")


            fmtStr = "{:.4f}"

            statdf = DataFrame()
            statdf.setTitle("<h3>Alignment Statistics</h3>")
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

            statdf.addRow(DataRow.fromDict({
                'Name': 'Avg GC content of Reference Sequence',
                'Value': fmtStr.format(gcContentRef/refLength)
            }))


            args.pltcfg.makeTable(statdf)



            allSubst = sum([ntSubstitutions[x] for x in ntSubstitutions])

            substDF = DataFrame()
            substDF.setTitle("<h3>Substitution Statistics (read -> ref)</h3>")
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
                kmerDF.setTitle("<h3>k-mer histogram: perfectly aligned (k={})</h3>".format(args.perfectk))
                args.pltcfg.makeTable(kmerDF)


                for evtype in kmerCounter:

                    obsInfo = {
                        'RUNID': fileName + "_" + str(evtype),
                        'USER_RUN_NAME': evtype,
                        'KMERCOUNTS': kmerCounter[evtype]
                    }

                    kmerDF = KmerHistogram.dfSummary(obsInfo, args.mc)
                    kmerDF.setTitle("<h3>k-mer histogram: sequence before CIGAR " + evtype + " (k={})</h3>".format(args.errork))
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
                kmerDF.setTitle("<h3>k-mer histogram: perfectly aligned (k={}, {})</h3>".format(args.perfectk, str(readType)))
                args.pltcfg.makeTable(kmerDF)

                for rtype in kmersByReadType[readType]:

                    print(rtype)

                    obsInfo = {
                        'RUNID': fileName + "_" + str(readType) + "_" + str(rtype),
                        'USER_RUN_NAME': rtype,
                        'KMERCOUNTS': kmersByReadType[readType][rtype]
                    }
                    kmerDF = KmerHistogram.dfSummary(obsInfo, args.mc)
                    kmerDF.setTitle("<h3>k-mer histogram: sequence before CIGAR " + rtype + " (k={}, {})</h3>".format(args.errork, str(readType)))
                    args.pltcfg.makeTable(kmerDF)

                overallIdentity = (sum(totalCounter[readType]["="]) / totalReadBases[readType]) * 100
                alignedIdentity = (sum(totalCounter[readType]["="]) / (sum(totalCounter[readType]['=']) + sum(totalCounter[readType]['X']))) * 100

                alignedIdentity = (sum(totalCounter[readType]["="]) / sum(totalCounter[readType]['M'])) * 100
                insertedPerAligned = (sum(totalCounter[readType]["I"]) / sum(totalCounter[readType]['M'])) * 100
                deletedPerAligned = (sum(totalCounter[readType]["D"]) / sum(totalCounter[readType]['M'])) * 100
                substitutedPerAligned = (sum(totalCounter[readType]["X"]) / sum(totalCounter[readType]['M'])) * 100

                fmtStr = "{:.4f}"

                statdf = DataFrame()
                statdf.setTitle("<h3>Alignment Statistics ({})</h3>".format(readType))
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
import pysam
import numpy as np
import HTSeq
import argparse
import time

import pickle
import sys
from collections import defaultdict
from porestat.utils.ArgParseExt import FileStubType, FolderType
from porestat.utils import eprint

from porestat.tools.kmer_coverage import KmerHistogram
from .ParallelAlignmentPSTReportableInterface import ParallelAlignmentPSTReportableInterface
from ..hdf5tool.Fast5File import Fast5TYPEAction
from ..plots.plotconfig import PlotConfig, PlotSaveTYPE
from ..hdf5tool import Fast5TYPE
from ..plots.poreplot import PorePlot

from ..utils.Files import fileExists, makePath
from ..utils.DataFrame import DataFrame, DataRow

from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException
from ..utils.Utils import mergeDicts
from ..utils.Stats import calcN50

from collections import Counter

import ctypes
from collections import Counter

import os

baseFolder = str(os.path.dirname(os.path.realpath(__file__)))

lib = ctypes.cdll.LoadLibrary(baseFolder+'/../../clib/lib/libfoo.so')


class COUNTER_PAIR(ctypes.Structure):
    _fields_ = [
        ("key", ctypes.c_char_p),
        ("value", ctypes.c_uint32)
    ]

class SeqCoverage:

    def __init__(self, seqid, start, end):

        self.seqid = seqid
        self.start = start
        self.end = end

class SEQ_COVERAGE(ctypes.Structure):
    _fields_ = [
        ("seqid", ctypes.c_uint32),
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

        print("WARNING: this is not designed to load multiple fasta!")

        print("Loading Fasta File", fastaPath)
        lib.AlignmentStatistics_load_fasta(self.obj, charStr)

    def processFiles(self, vals):

        sam2header = {}

        for x in vals:
            reader = pysam.Samfile(x.name, "rb", check_sq=False)
            samHeader = reader.header

            sam2header[x.name] = samHeader.to_dict()
            print(samHeader.to_dict())


        vals = [x.name.encode('utf-8') for x in vals]
        c_array = (ctypes.c_char_p * len(vals))(*vals)

        charPArray = ctypes.cast(c_array, ctypes.POINTER(ctypes.c_char_p))

        retSize = lib.AlignmentStatistics_process(self.obj, ctypes.c_int(len(vals)), charPArray)

        retVecSize = lib.AlignmentStatistics_ReadStats_size(self.obj)

        retVec = lib.AlignmentStatistics_ReadStats(self.obj)

        s = ctypes.cast(retVec, ctypes.POINTER(READ_ALIGNMENT * retVecSize))

        allElems = [x for x in s.contents]

        print("Processing C++ Results")

        readResults = {}

        for i in range(0, retVecSize):

            elem = allElems[i]

            """

            CIGAR2LEN
            """
            cigars = [x for x in elem.CIGARS[0:elem.CIGAR_COUNT].decode()]
            vecLengths = [x for x in ctypes.cast(elem.CIGAR_VEC_LENGTHS,
                                                 ctypes.POINTER(ctypes.c_uint32 * elem.CIGAR_COUNT)).contents]
            cigarVecs = [x for x in
                         ctypes.cast(elem.CIGAR_LENGTHS, ctypes.POINTER(ctypes.c_void_p * elem.CIGAR_COUNT)).contents]

            cigar2len = {}
            for c in range(0, len(cigars)):
                cigarLen = [x for x in
                            ctypes.cast(cigarVecs[c], ctypes.POINTER(ctypes.c_uint32 * vecLengths[c])).contents]

                cigar2len[cigars[c]] = cigarLen


            #if 'M' in cigar2len and 'Z' in cigar2len and 'E' in cigar2len:
            #    print(elem.pReadID.decode())
            #    print(sum(cigar2len['M']))
            #    print(sum(cigar2len['E']))
            #    print(sum(cigar2len['Z']))
            #       assert( sum(cigar2len['M']) == sum(cigar2len['E']) + sum(cigar2len['Z']) )


            """
 
            PERF KMERS
            """
            perfKmers = [x for x in
                         ctypes.cast(elem.PERF_KMERS, ctypes.POINTER(COUNTER_PAIR * elem.PERF_KMER_COUNT)).contents]

            perfKmerCounter = Counter()
            for x in perfKmers:
                perfKmerCounter[x.key.decode()] = x.value


            """

            CIGAR2LEN
            """
            mmLengths = [x for x in ctypes.cast(elem.MM2KMER_VEC_LENGTHS,
                                                ctypes.POINTER(ctypes.c_uint32 * elem.CIGAR_COUNT)).contents]
            mm2kmers = [x for x in
                        ctypes.cast(elem.MM2KMER_LENGTHS, ctypes.POINTER(ctypes.c_void_p * elem.CIGAR_COUNT)).contents]

            mm2kmer = defaultdict(lambda: Counter())
            for c in range(0, len(cigars)):
                mmKmers = [x.decode() for x in
                           ctypes.cast(mm2kmers[c], ctypes.POINTER(ctypes.c_char_p * mmLengths[c])).contents]

                for kmer in mmKmers:
                    mm2kmer[cigars[c]][kmer] += 1


            """

            COVERAGES
            """
            covIntervals = [x for x in ctypes.cast(elem.COVERAGES, ctypes.POINTER(SEQ_COVERAGE * elem.COV_COUNT)).contents]
            coverageIntervals = []

            for cov in covIntervals:
                coverageIntervals.append(SeqCoverage(cov.seqid, cov.start, cov.end))

            # TODO add this to global array

            """

            NT SUBSTITUTION
            """
            ntSubst = [x for x in
                       ctypes.cast(elem.NT_SUBST, ctypes.POINTER(NT_SUBSTITUTION * elem.NT_SUBST_COUNT)).contents]

            substCounter = Counter()
            for x in ntSubst:
                #print(x)
                #print(x.src, chr(x.src))
                #print(x.tgt, chr(x.tgt))
                #print(x.count)
                substCounter[(x.src.decode(), x.tgt.decode())] = x.count


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
            readStat.coverages = coverageIntervals
            readStat.mm2kmer = mm2kmer

            readResults[readStat.read_id] = readStat

        print("Finished C++ Results", len(readResults), len(sam2header))

        return readResults, sam2header



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

        #parser.add_argument('-q', '--read-type', nargs='+', dest='read_type', action=Fast5TYPEAction, default=None)
        parser.add_argument('-v', '--violin', dest='violin', action='store_true', default=False)

        parser.add_argument('-o', '--output', type=FolderType('w'), help='output folder for report', required=True)
        parser.add_argument('-n', '--output-name', type=str, help='output name', required=False)

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

        self.genomicArray = HTSeq.GenomicArray( [], stranded=False )
        self.sam2header = {}

        self.args.pltcfg.setOutputType(PlotSaveTYPE.HTML_STRING)



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

    def getRead2Type(self):
        
        r2t = {}

        if self.hasArgument('read_info', self.args) and self.args.readInfo:
            idIdx = self.args.readInfo.getColumnIndex("READ_ID")
            typeIdx = self.args.readInfo.getColumnIndex("TYPE")

            for row in self.args.readInfo.data:
                r2t[ row[idIdx] ] = row[typeIdx]


        return r2t


    def getReadInfo(self, readName):

        if self.hasArgument('read_info', self.args) and self.args.readInfo:

            readInfo = self.args.readInfo.findRow('READ_ID', readName)
            if readInfo != None:
                readType = Fast5TYPE[ readInfo[ "TYPE" ] ]
                readID = readInfo[ "READ_ID" ]

            return readType,readID

        return "Unknown", readName


    def addCoverage(self, chrom, start, end):

        if not chrom in self.genomicArray.chrom_vectors:
            self.genomicArray.add_chrom(chrom)

        self.genomicArray[HTSeq.GenomicInterval( chrom, start, end, "." )] += 1

    def joinParallel(self, existResult, newResult, oEnvironment):

        print("Joining Results")
        
        if existResult is None:
            existResult = {"align": {}, "sam2header": {}}

        alignRes = newResult[0]
        sam2hRes = newResult[1]

        existResult["align"].update(alignRes)
        existResult["sam2header"].update(sam2hRes)

        return existResult


    def execParallel(self, data, environment):
        """

        :param data: input files
        :param environment: usually NONE
        :return: None
        """
        statMaker = CAlignmentStatistics()

        for fastaFile in self.args.fasta_files:
            statMaker.loadFasta(fastaFile.name)


        alignResults = {}
        sam2HeaderResults = {}

        for samFile in data:
            print("Processing samFile", samFile.name)
            readStats, sam2header = statMaker.processFiles([samFile])
            print("Received data", len(readStats))

            for x in sam2header:
                sam2HeaderResults[x] = sam2header[x]

            print("Applying Read Info")
            readId2Type = self.getRead2Type()
            for rid,readID in enumerate(readStats):
                if rid % 100000 == 0:
                    print("Processing read", rid)
                readStats[readID].read_type = readId2Type.get(readID, "Unknown")

            print("got", len(readStats), "readstats")
            print("Finished Processing", samFile.name)
            alignResults[samFile.name] = readStats


        return {"align": alignResults, "sam2header": sam2HeaderResults}


    def exec(self):

        iStart = time.time()
        inputs = self.prepareInputs(self.args)
        iEnd = time.time()
        eprint("Preparing Inputs: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd - iStart))))

        print(inputs)
        environment = self.prepareEnvironment(self.args)
        llResults = self.execParallel(inputs, environment)
        iEnd = time.time()
        eprint("Exec Parallel: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd - iStart))))

        self.makeResults(llResults, None, self.args)
        iEnd = time.time()
        eprint("Results Time: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd-iStart))))


    def getSeqNameFromId(self, samFileName, seqid):
        return self.sam2header[samFileName]["SQ"][seqid]["SN"]


    @classmethod
    def mean(cls, elems):

        return sum(elems) / len(elems)

    def makeResults(self, alignSam2HRes, oEnvironment, args):
        """
        :param allFilesResult:
        :param oEnvironment:
        :param args:
        :return:
        """

        allFilesResult = alignSam2HRes["align"]
        self.sam2header = alignSam2HRes["sam2header"]


        print("Preparing results")

        args.output = makePath(args.output)

        print("Output folder: " + str(args.output))
        print("Output name:   " + str(args.output_name))

        if args.output_name == None:
            args.output_name = 'report'
            print("Changed Output name:   " + str(args.output_name))

        """
        
        FASTA REFERENCE INFOS
        
        """

        gcContentRef = 0
        refLength = 0

        fastaSeqs = {}
        for fastaFile in self.args.fasta_files:

            for seq in HTSeq.FastaReader(fastaFile):
                fastaSeqs[seq.name] = str(seq.seq)

                strSeq = str(seq.seq)
                gcContentRef += self.calc_gc(strSeq) * len(strSeq)
                refLength += len(strSeq)



        """
        
        GTF/FEATURE INFOS
        
        """
        allFeatureArrays = defaultdict(lambda: HTSeq.GenomicArray([], stranded=False))

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

        readCoverage = HTSeq.GenomicArray([], stranded=False)

        for fileName in allFilesResult:



            args.pltcfg.addHTMLPlot("<h2>" + str(fileName) + "</h2>\n")
            args.pltcfg.addHTMLPlot("<h3>All Reads Overview</h3>\n")

            allReads = allFilesResult[fileName]

            alignedIDs = set()
            unalignedIDs = set()

            allAlignedReads = {}
            for x in allReads:

                elem = allReads[x]

                if elem.aligned == True:

                    allAlignedReads[x] = elem
                    alignedIDs.add(elem.read_id)

                    for cov in elem.coverages:

                        covChr = self.getSeqNameFromId(fileName, cov.seqid)

                        if not covChr in readCoverage.chrom_vectors:
                            readCoverage.add_chrom(covChr)

                        #print(cov.seq, cov.start, cov.end)
                        readCoverage[HTSeq.GenomicInterval(covChr, cov.start, cov.end+1, ".")] += 1


                else:
                    unalignedIDs.add(elem.read_id)



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

            print("Reads by status", plotDataDict)

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

            #kmerCounter = defaultdict(lambda: Counter())

            totalCounter = defaultdict(list)
            perfectCounter = Counter()

            allReadInfos = defaultdict(list)
            totalReadBases = 0

            ntSubstitutions = Counter()

            alignedIdx = 0
            
            allReadInfos['READ_IDENTITY'] = [(None, None)] * len(alignedIDs)
            allReadInfos['ALIGNMENT_IDENTITY'] = [(None, None)] * len(alignedIDs)
            allReadInfos['MATCHED_LONGEST_LENGTH'] = [(None, None)] * len(alignedIDs)
            allReadInfos['READ_MATCHED_GC_CONTENT'] = [(None, None)] * len(alignedIDs)
            allReadInfos['REF_MATCHED_GC_CONTENT'] = [(None, None)] * len(alignedIDs)
            allReadInfos['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH'] = [(None, None)] * len(alignedIDs)
            allReadInfos['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH'] = [(None, None)] * len(alignedIDs)

            allReadInfos['READ_ALIGN_QUAL_TO_READ_LENGTH'] = [(None, None)] * len(alignedIDs)
            allReadInfos['SEQ_QUAL_MIN_TO_READ_LENGTH'] = [(None, None)] * len(alignedIDs)
            allReadInfos['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN'] = [(None, None)] * len(alignedIDs)
            allReadInfos['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN'] = [(None, None)] * len(alignedIDs)
            
            for ridx, x in enumerate(allReads):

                if ridx % 100000 == 0:
                    print("Processing Read", ridx)

                elem = allReads[x]

                if not elem.aligned:
                    continue


                """
                {
                    'READ_IDENTITY': (readIdentityCount / len(readAlignment.read), len(readAlignment.read)),
                    'ALIGNMENT_IDENTITY': (alignmentIdentityCount / alignmentLength, alignmentLength),
                    'MATCHED_GC_CONTENT': matchedGCContents,
                    'MATCHED_LONGEST_LENGTH': [matchedLongestLength]
                }
                """

                allReadInfos['READ_IDENTITY'][alignedIdx] = (elem.read_identity, elem.read_length)
                allReadInfos['ALIGNMENT_IDENTITY'][alignedIdx] = (elem.ref_identity, elem.read_length)
                allReadInfos['MATCHED_LONGEST_LENGTH'][alignedIdx] = elem.longest_matched_stretch

                allReadInfos['READ_MATCHED_GC_CONTENT'][alignedIdx] = (elem.read_gc, sum(elem.cigar_lengths['M'])) 
                allReadInfos['REF_MATCHED_GC_CONTENT'][alignedIdx] = (elem.ref_gc, sum(elem.cigar_lengths['M'])) 


                allReadInfos['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH'][alignedIdx] = ( elem.read_gc  , elem.read_length)  
                allReadInfos['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH'][alignedIdx] = (  elem.ref_gc, elem.ref_length)     

                totalReadBases += elem.read_length

                allReadInfos['READ_ALIGN_QUAL_TO_READ_LENGTH'][alignedIdx] = (elem.align_quality, elem.read_length)   
                allReadInfos['SEQ_QUAL_MIN_TO_READ_LENGTH'][alignedIdx] = (elem.seqqual_min, elem.read_length)     
                allReadInfos['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN'][alignedIdx] = (elem.align_quality, elem.seqqual_min)   
                allReadInfos['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN'][alignedIdx] = (elem.align_quality, elem.seqqual_median)

                perfectCounter.update(elem.perfect_kmers)

                for cl in elem.cigar_lengths:
                    totalCounter[cl] += elem.cigar_lengths[cl]

                #for x in elem.perfect_kmers:
                #    kmerCounter[x] += elem.perfect_kmers[x]

                ntSubstitutions.update(elem.nt_substitutions)
                #for mm in elem.nt_substitutions:
                #    ntSubstitutions[mm] += elem.nt_substitutions[mm]

                alignedIdx += 1

            # I = Insertion
            # D = Deletion
            # M = Match (either exact or not)
            # E = Exact Match
            # Z = Inexact Match

            overallIdentity = (sum(totalCounter["E"]) / totalReadBases) * 100
            alignedMatch = (sum(totalCounter["M"]) / (totalReadBases)) * 100
            matchedIdentity = (sum(totalCounter["E"]) / sum(totalCounter['M'])) * 100

            insertedPerAligned = (sum(totalCounter["I"]) / sum(totalCounter['M'])) * 100
            deletedPerAligned = (sum(totalCounter["D"]) / sum(totalCounter['M'])) * 100
            substitutedPerAligned = (sum(totalCounter["Z"]) / sum(totalCounter['M'])) * 100




            args.pltcfg.addHTMLPlot("<h2>Feature-based Coverage</h2>\n")

            if args.gtf and len(allFeatureArrays) > 0:


                fmtStr = "{:.4f}"

                featureDF = DataFrame()
                featureDF.setTitle("<h3>Alignment Feature Statistics</h3>")
                featureDF.addColumns(['Name', 'Value'])



                for ftype in allFeatureArrays:

                    featureCovered = 0
                    readsCovered = 0

                    missingChromsInReadCoverage = set()


                    fga = allFeatureArrays[ftype]

                    for featurePart in fga.steps():

                        partCovered = featurePart[0].length
                        partCount = featurePart[1]

                        if partCount > 0:

                            featureCovered += partCovered

                            #for this part fetch steps in read coverage
                            if not featurePart[0].chrom in readCoverage.chrom_vectors:
                                missingChromsInReadCoverage.add(featurePart[0].chrom)
                                continue

                            for readPart in readCoverage[featurePart[0]].steps():

                                readCovered = readPart[0].length
                                readCount = readPart[1]

                                if readCount > 0:
                                    readsCovered += readCovered

                    print("Missing Chromosomes in ReadCoverage:",  missingChromsInReadCoverage)

                    featureDF.addRow(DataRow.fromDict({
                        'Name': 'Featuretype Coverage ({})'.format(ftype),
                        'Value': fmtStr.format(featureCovered/refLength)
                    }))

                    featureDF.addRow(DataRow.fromDict({
                        'Name': 'Featuretype Covered By Reads ({})'.format(ftype),
                        'Value': fmtStr.format(readsCovered/featureCovered)
                    }))

                args.pltcfg.makeTable(featureDF)

            else:

                args.pltcfg.addHTMLPlot("<p>No feature annotation supplied.</p>\n")


            args.pltcfg.addHTMLPlot("<h3>Alignment Statistics</h3>\n")


            fmtStr = "{:.4f}"


            statdf = DataFrame()
            statdf.setTitle("<h3>CIGAR counts</h3>")
            statdf.addColumns(['Name', 'Value'])

            statdf.addRow(DataRow.fromDict({
                'Name': "Total",
                'Value': fmtStr.format(totalReadBases)
            }))

            for cigarElem in totalCounter:
                statdf.addRow(DataRow.fromDict({
                    'Name': cigarElem,
                    'Value': fmtStr.format(sum(totalCounter[cigarElem]))
                }))

            args.pltcfg.makeTable(statdf)
            print(statdf)

            statdf = DataFrame()
            statdf.setTitle("<h3>Alignment Statistics</h3>")
            statdf.addColumns(['Name', 'Value'])

            statdf.addRow(DataRow.fromDict({
                'Name': 'Overall Identity (Exact/All)',
                'Value': fmtStr.format(overallIdentity)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Matching Bases (M/All)',
                'Value': fmtStr.format(alignedMatch)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Matched Identity/Identical bases per 100 aligned bases (E/M)',
                'Value': fmtStr.format(matchedIdentity)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Inserted bases per 100 aligned bases (I/M)',
                'Value': fmtStr.format(insertedPerAligned)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Deleted bases per 100 aligned bases (D/M)',
                'Value': fmtStr.format(deletedPerAligned)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Substituted bases per 100 aligned bases (Z/M)',
                'Value': fmtStr.format(substitutedPerAligned)
            }))

            addedGC = sum([x[0]*x[1] for x in allReadInfos['READ_MATCHED_GC_CONTENT']])
            addedLength = sum([x[1] for x in allReadInfos['READ_MATCHED_GC_CONTENT']])

            statdf.addRow(DataRow.fromDict({
                'Name': 'Avg GC content of Matched/Aligned Read Sequence',
                'Value': fmtStr.format(100.0*addedGC/addedLength)
            }))

            statdf.addRow(DataRow.fromDict({
                'Name': 'Avg GC content of Reference Sequence',
                'Value': fmtStr.format(100.0*gcContentRef/refLength)
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



            if self.hasArgument('read_info', args) and args.readInfo:

                obsInfo = {
                    'RUNID': fileName + "_perfect",
                    'USER_RUN_NAME': "perfect",
                    'KMERCOUNTS': perfectCounter
                }

                kmerDF = KmerHistogram.dfSummary(obsInfo, args.mc)
                kmerDF.setTitle("<h3>k-mer histogram: perfectly aligned (k={})</h3>".format(args.perfectk))
                args.pltcfg.makeTable(kmerDF)

                """
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




                """
                READ IDENTITY, ALIGNMENT IDENTITY, MATCHED GC CONTENT, MAX MATCHED PART
                """

                xReadID = [x[0] for x in allReadInfos['READ_IDENTITY']]
                yReadID = [x[1] for x in allReadInfos['READ_IDENTITY']]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Identity vs Read Length", "Read Identity", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos['ALIGNMENT_IDENTITY']]
                yReadID = [x[1] for x in allReadInfos['ALIGNMENT_IDENTITY']]
                print("Mean Align Identity", np.mean(xReadID), np.median(xReadID))
                print("Mean Read Length", np.mean(yReadID), np.median(yReadID))
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Identity vs Read Length", "Alignment Identity", "Read Length [bp]", pltcfg=args.pltcfg)

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

            allReadTypes = set()
            alignedByReadType = Counter()
            for ridx, x in enumerate(allReads):

                if ridx % 100000 == 0:
                    print("Precalculating by Read Type; Processing Read", ridx)

                elem = allReads[x]
                readType = str(elem.read_type)
                allReadTypes.add(readType)

                if elem.aligned:
                    alignedByReadType[readType] += 1


            for readType in allReadTypes:

                allReadInfos[readType]['READ_IDENTITY'] = [(None, None)] * (alignedByReadType[readType])
                allReadInfos[readType]['ALIGNMENT_IDENTITY'] = [(None, None)] * (alignedByReadType[readType])
                allReadInfos[readType]['MATCHED_LONGEST_LENGTH'] = [(None, None)] * (alignedByReadType[readType])
                allReadInfos[readType]['READ_MATCHED_GC_CONTENT'] = [(None, None)] * (alignedByReadType[readType])
                allReadInfos[readType]['REF_MATCHED_GC_CONTENT'] = [(None, None)] * (alignedByReadType[readType])
                allReadInfos[readType]['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH'] = [(None, None)] * (alignedByReadType[readType])
                allReadInfos[readType]['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH'] = [(None, None)] * (alignedByReadType[readType])

                allReadInfos[readType]['READ_ALIGN_QUAL_TO_READ_LENGTH'] = [(None, None)] * (alignedByReadType[readType])
                allReadInfos[readType]['SEQ_QUAL_MIN_TO_READ_LENGTH'] = [(None, None)] * (alignedByReadType[readType])
                allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN'] = [(None, None)] * (alignedByReadType[readType])
                allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN'] = [(None, None)] * (alignedByReadType[readType])


            totalReadBases = Counter()

            alignedIdx = Counter()

            for ridx, x in enumerate(allReads):

                if ridx % 100000 == 0:
                    print("Processing Read", ridx)

                elem = allReads[x]

                if not elem.aligned:
                    continue

                readType = str(elem.read_type)

                #allReadInfos[readType] = mergeDicts(allReadInfos[readType], allReadInfo)

                allReadInfos[readType]['READ_IDENTITY'][alignedIdx[readType]] = (elem.read_identity, elem.read_length)
                allReadInfos[readType]['ALIGNMENT_IDENTITY'][alignedIdx[readType]] = (elem.ref_identity, elem.read_length)

                totalReadBases[readType] += elem.read_length

                allReadInfos[readType]['READ_MATCHED_GC_CONTENT'][alignedIdx[readType]] = (elem.read_gc, sum(elem.cigar_lengths['M']))
                allReadInfos[readType]['REF_MATCHED_GC_CONTENT'][alignedIdx[readType]] = (elem.ref_gc, sum(elem.cigar_lengths['M']))

                allReadInfos[readType]['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH'][alignedIdx[readType]] = ( elem.read_gc  , elem.read_length)
                allReadInfos[readType]['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH'][alignedIdx[readType]] = ( elem.ref_gc  , elem.ref_length)

                allReadInfos[readType]['READ_ALIGN_QUAL_TO_READ_LENGTH'][alignedIdx[readType]] = (elem.align_quality, elem.read_length)
                allReadInfos[readType]['SEQ_QUAL_MIN_TO_READ_LENGTH'][alignedIdx[readType]] = (elem.seqqual_min, elem.read_length)
                allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN'][alignedIdx[readType]] = (elem.align_quality, elem.seqqual_min)

                allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN'][alignedIdx[readType]] = (elem.align_quality, elem.seqqual_median)

                alignedIdx[readType] += 1

                perfectCounter[readType].update(elem.perfect_kmers)

                for cl in elem.cigar_lengths:
                    totalCounter[readType][cl] += elem.cigar_lengths[cl]

                kmersByReadType[readType].update(elem.mm2kmer)

                #for x in elem.nt_substitutions:
                #    ntSubstitutions[read_type][x] += elem.nt_substitutions[x]

            print(sum([allReads[x].read_length for x in allReads if
                 allReads[x].aligned and allReads[x].read_type == Fast5TYPE.BASECALL_2D]))

            for readType in totalReadBases:
                print(readType, totalReadBases[x])

            for readType in totalCounter:

                print(readType)
                print(self.hasArgument('violin', args),args.violin)

                if totalReadBases[readType] == 0:
                    print("Skipping read type ", readType, " for having no bases sequenced")
                    continue

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

                overallIdentity = (sum(totalCounter[readType]["E"]) / totalReadBases[readType]) * 100
                alignedIdentity = (sum(totalCounter[readType]["E"]) / (sum(totalCounter[readType]['E']) + sum(totalCounter[readType]['Z']))) * 100

                alignedIdentity = (sum(totalCounter[readType]["E"]) / sum(totalCounter[readType]['M'])) * 100
                insertedPerAligned = (sum(totalCounter[readType]["I"]) / sum(totalCounter[readType]['M'])) * 100
                deletedPerAligned = (sum(totalCounter[readType]["D"]) / sum(totalCounter[readType]['M'])) * 100
                substitutedPerAligned = (sum(totalCounter[readType]["Z"]) / sum(totalCounter[readType]['M'])) * 100

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

                xReadID = [x[0] for x in allReadInfos[readType]['READ_IDENTITY'] if not None in x]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_IDENTITY'] if not None in x]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Identity vs Read Length ({})".format(str(readType)), "Read Identity", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['ALIGNMENT_IDENTITY'] if not None in x]
                yReadID = [x[1] for x in allReadInfos[readType]['ALIGNMENT_IDENTITY'] if not None in x]
                print("Mean Align Identity", np.mean(xReadID), np.median(xReadID))
                print("Mean Read Length", np.mean(yReadID), np.median(yReadID))
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Identity vs Read Length  ({})".format(str(readType)), "Alignment Identity", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['READ_MATCHED_GC_CONTENT'] if not None in x]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_MATCHED_GC_CONTENT'] if not None in x]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Matched GC Content vs Aligned Seq Length ({})".format(str(readType)), "GC Content", "Aligned Sequence Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['REF_MATCHED_GC_CONTENT'] if not None in x]
                yReadID = [x[1] for x in allReadInfos[readType]['REF_MATCHED_GC_CONTENT'] if not None in x]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Reference Matched GC Content vs Aligned Seq Length ({})".format(str(readType)), "GC Content", "Aligned Sequence Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH'] if not None in x]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_MATCHED_GC_CONTENT_TO_READ_LENGTH'] if not None in x]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Matched GC Content vs Read Length ({})".format(str(readType)), "GC Content", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH'] if not None in x]
                yReadID = [x[1] for x in allReadInfos[readType]['REF_MATCHED_GC_CONTENT_TO_READ_LENGTH'] if not None in x]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Reference Matched GC Content vs Read Length ({})".format(str(readType)), "GC Content", "Read Length [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_READ_LENGTH'] if not None in x]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_READ_LENGTH'] if not None in x]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Read Align Quality vs Read Length ({})".format(str(readType)), "Read Align Quality", "Read Lenght [bp]", pltcfg=args.pltcfg)

                xReadID = [x[0] for x in allReadInfos[readType]['SEQ_QUAL_MIN_TO_READ_LENGTH'] if not None in x]
                yReadID = [x[1] for x in allReadInfos[readType]['SEQ_QUAL_MIN_TO_READ_LENGTH'] if not None in x]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Minimum Sequence Quality vs Read Length (min value of read, {})".format(str(readType)), "Minimum Sequence Quality", "Read Length [bp]", pltcfg=args.pltcfg)


                xReadID = [x[0] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN'] if not None in x]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MIN'] if not None in x]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Quality vs Minimum Sequence Quality ({})".format(str(readType)), "Alignment Quality", "Sequence Quality (min value of read)", pltcfg=args.pltcfg)


                xReadID = [x[0] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN'] if not None in x]
                yReadID = [x[1] for x in allReadInfos[readType]['READ_ALIGN_QUAL_TO_SEQ_QUAL_MEAN'] if not None in x]
                PorePlot.plot_scatter_densities(xReadID, yReadID, "Alignment Quality vs Mean Sequence Quality", "Alignment Quality", "Sequence Quality (mean value of read)", pltcfg=args.pltcfg)



                if self.hasArgument('violin', args) and args.violin:
                    PorePlot.plotViolin(totalCounter[readType], sorted([x for x in totalCounter[readType]]), "CIGARs: " + str(readType),
                                        pltcfg=args.pltcfg, shareX = False, shareY=False)
                else:
                    PorePlot.plotBoxplot(totalCounter[readType], sorted([x for x in totalCounter[readType]]), "CIGARs: " + str(readType),
                                         pltcfg=args.pltcfg)



        args.pltcfg.prepareHTMLOutput(args.output, args.output_name + ".html", relativeImport=True)

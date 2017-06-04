import argparse
import HTSeq


from ..utils.DataFrame import DataFrame
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Parallel import Parallel as ll
from ..utils.Utils import mergeDicts, fileExists
from ..utils.Stats import calcN50


class ReadCountAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(ReadCountAnalysisFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('counts', help='expls help')
        parser.add_argument('-s', '--sam', nargs='+', type=str, required=True, help='alignment files')
        parser.add_argument('-g', '--gff', type=argparse.FileType("r"), required=True, help='gene annotation')
        parser.add_argument('-r', '--read_info', nargs='+', type=str, help='read summary file', required=False)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return ReadCountAnalysis(simArgs)



class ReadCountAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(ReadCountAnalysis, self).__init__(args)

        self.readInfo = None
        self.genomeAnnotation = None

    def _makePropDict(self):

        props = {}

        props['COUNTS'] = {}
        props['COUNTS']['ALL'] = 0

        return props

    def readGenomeAnnotation(self, args):

        self.genomeAnnotation = HTSeq.GFF_Reader( args.gff, end_included=True )

    def readReadInfo(self, args):

        self.readInfo = {}

        if args.read_info == None:
            return

        if not fileExists(args.read_info):
            PSToolException("Read info file does not exist: " + args.read_info)

        allReadData = DataFrame.parseFromFile(args.read_info, cDelim='\t')

        readNameIdx = allReadData.getColumnIndex("READ_NAME")
        readTypeIdx = allReadData.getColumnIndex("TYPE")
        readRunIdx = allReadData.getColumnIndex("RUN_NAME")
        readLengthIdx = allReadData.getColumnIndex("READ_LENGTH")

        for elem in allReadData.vElements:

            readName = elem[readNameIdx]
            readType = elem[readTypeIdx]
            readRun  = elem[readRunIdx]
            readLength = int(elem[readLengthIdx])

            self.readInfo[readName] = (readType, readRun, readLength)

    def prepareInputs(self, args):
        self.readReadInfo(args)

    def execParallel(self, data, environment):

        cvg = HTSeq.GenomicArray("auto", stranded=False, typecode='i')
        alignment_file = HTSeq.SAM_Reader( environment.sam )

        foundReadsAligned = set()
        foundReadsNotAligned = set()

        for alngt in alignment_file:
            if alngt.aligned:

                foundReadsAligned.add( alngt.read.name )

                cigars = alngt.cigar

                for cigar in cigars:

                    if cigar.type == 'M':
                        cvg[cigar.ref_iv] += 1

            else:
                foundReadsNotAligned.add( alngt.read.name )

        featureLengths = {}
        featureCoverage = {}
        for feature in self.genomeAnnotation:

            if feature.type == "exon":

                if not feature.name in featureLengths:
                    featureLengths[feature.name] = []

                featureLengths[feature.name].append(feature.iv.end - feature.iv.start)

                if not feature.name in featureCoverage:
                    featureCoverage[feature.name] = []

                featureCoverage[feature.name].append(cvg[feature.iv])

        allKeys = []

        for x in featureCoverage:
            allKeys.append(x)

        coverages = Counter()
        covlengths = Counter()

        summedAverage = Counter()

        for x in sorted(allKeys):
            featureCount = 0

            for y in featureCoverage[x]:

                for iv, value in y.steps():
                    featureCount += value * (iv.end - iv.start)

            # print str(x) + " " + str(featureCount) + " " + str(sum(featureLengths[x]))

            featureLength = sum(featureLengths[x])
            cov = float(featureCount) / float(featureLength)
            print(str(x) + " " + str(cov))

            if x.split("_")[1].startswith("r"):
                coverages["ribo"] += featureCount
                covlengths["ribo"] += featureLength

                summedAverage["ribo"] += cov
            else:
                coverages["nonribo"] += featureCount
                covlengths["nonribo"] += featureLength

                summedAverage["nonribo"] += cov

        for x in coverages:
            print(str(x) + " " + str(float(coverages[x]) / float(covlengths[x])))
            print(str(x) + " " + str(summedAverage))


    def joinParallel(self, existResult, newResult, oEnvironment):

        return None


    def makeResults(self, parallelResult, oEnvironment, args):

        pass



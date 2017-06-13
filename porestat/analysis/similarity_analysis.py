import argparse
import HTSeq
from porestat.plots.poreplot import PorePlot

from porestat.plots.plotconfig import PlotConfig
from ..utils.DataFrame import DataFrame
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Files import fileExists
from scipy import stats
import sys
import numpy as np


class SimilarityAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(SimilarityAnalysisFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('counts', help='expls help')
        parser.add_argument('-s', '--sam', nargs='+', type=str, required=True, help='alignment files')
        parser.add_argument('-g', '--gff', type=argparse.FileType("r"), required=True, help='gene annotation')
        parser.add_argument('-r', '--read-info', nargs='+', type=str, help='read summary file', required=False)
        parser.add_argument('--no-plot', dest='plot', action='store_true', default=False)

        def fileOpener( filename ):

            open(filename, 'w').close()

            return filename

        parser.add_argument('-o', '--output', type=fileOpener, help='output location, default: std out', default=sys.stdout)

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return SimilarityAnalysis(simArgs)



class SimilarityAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(SimilarityAnalysis, self).__init__(args)

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

        for readInfoFile in args.read_info:

            if not fileExists(readInfoFile):
                raise PSToolException("Read info file does not exist: " + str(readInfoFile))

        self.readInfo = {}

        for readInfoFile in args.read_info:
            allReadData = DataFrame.parseFromFile(readInfoFile, cDelim='\t')

            readNameIdx = allReadData.getColumnIndex("READ_NAME")
            readTypeIdx = allReadData.getColumnIndex("TYPE")
            readRunIdx = allReadData.getColumnIndex("USER_RUN_NAME")
            readLengthIdx = allReadData.getColumnIndex("READ_LENGTH")

            for elem in allReadData.vElements:
                readName = elem[readNameIdx]
                readType = elem[readTypeIdx]
                readRun = elem[readRunIdx]
                readLength = int(elem[readLengthIdx])

                self.readInfo[readName] = (readType, readRun, readLength)

    def prepareInputs(self, args):
        self.readReadInfo(args)
        self.readGenomeAnnotation(args)

        for x in args.sam:
            if not fileExists(x):
                PSToolException("sam file does not exist: " + str(x))

        self.writeLinesToOutput(args.output, "\t".join(['gene', 'coverage', 'rank']) + "\n", mode='w')


        return args.sam

    def execParallel(self, data, environment):

        cvg = HTSeq.GenomicArray("auto", stranded=False, typecode='i')
        alignment_file = HTSeq.SAM_Reader( data )

        foundReadsAligned = set()
        foundReadsNotAligned = set()

        for alngt in alignment_file:

            if not alngt.read.name in self.readInfo:
                self.readInfo[alngt.read.name] = (Fast5TYPE.UNKNOWN, None, len(alngt.read.seq))


            if alngt.aligned:
                foundReadsAligned.add( alngt.read.name )
                cigars = alngt.cigar
                for cigar in cigars:

                    if cigar.type == 'M':
                        cvg[cigar.ref_iv] += 1

            else:
                foundReadsNotAligned.add( alngt.read.name )

        print("Found aligned reads: " + str(len(foundReadsAligned)))
        print("Found unaligned reads: " + str(len(foundReadsNotAligned)))
        print("Total reads: " + str(len(foundReadsAligned) + len(foundReadsNotAligned)))

        totalBaseCount = sum([ self.readInfo[x][2] for x in self.readInfo  ])
        print("Total bases in files: " + str(totalBaseCount))

        featureLengths = {}
        featureCoverage = {}
        for feature in self.genomeAnnotation:

            if feature.type == "exon" or feature.type == "operon":

                featureLocusName = feature.attr['gene_id'] if 'gene_id' in feature.attr else feature.name

                if not featureLocusName in featureLengths:
                    featureLengths[ featureLocusName ] = []

                featureLengths[ featureLocusName ].append(feature.iv.end - feature.iv.start)

                if not feature.name in featureCoverage:
                    featureCoverage[ featureLocusName ] = []

                featureCoverage[ featureLocusName ].append(cvg[feature.iv])

        allKeys = []

        for x in featureCoverage:
            allKeys.append(x)

        coverages = Counter()
        covlengths = Counter()

        summedAverage = Counter()

        allCovs = []

        for x in sorted(allKeys):
            featureCount = 0

            for y in featureCoverage[x]:

                for iv, value in y.steps():
                    featureCount += value * (iv.end - iv.start)

            # print str(x) + " " + str(featureCount) + " " + str(sum(featureLengths[x]))

            featureLength = sum(featureLengths[x])
            cov = float(featureCount) / float(featureLength)
            #print(str(x) + " " + str(cov))

            if not x.startswith("operon"):
                allCovs.append( (x, cov) )

            xa = x.split("_")
            if len(xa) > 1 and xa[1].startswith("r"):
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

        ranked = stats.rankdata( [x[1] for x in allCovs] )

        allRankedCovs = []
        for i in range(0, len(allCovs)):
            allRankedCovs.append( (allCovs[i][0], allCovs[i][1], ranked[i]) )


        allLines = []
        for x in sorted(allRankedCovs, key=lambda x: x[2]):
            allLines.append( str(x[0]) + "\t" + str(x[1]) + "\t" + str(x[2]) + "\n" )

        self.writeLinesToOutput(environment.output, allLines)

        return allRankedCovs


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = []

        return existResult + newResult


    def makeResults(self, parallelResult, oEnvironment, args):

        if args.plot:

            plotData = {}
            plotData['coverage'] = [x[1] for x in parallelResult]

            PorePlot.plotCumHistogram( plotData, [x for x in plotData], 'Coverage Plot', bins=len(parallelResult), pltcfg=args.pltcfg )






import argparse
import HTSeq
from porestat.plots.poreplot import PorePlot

from porestat.plots.plotconfig import PlotConfig
from ..utils.DataFrame import DataFrame
from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter, defaultdict
from ..utils.Files import fileExists
from scipy import stats
import sys
import numpy as np


class ReadCountAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(ReadCountAnalysisFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-s', '--sam', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
        parser.add_argument('-g', '--gff', type=argparse.FileType("r"), required=True, help='gene annotation')
        parser.add_argument('-r', '--read-info', nargs='+', type=argparse.FileType('r'), help='read summary file', required=False)
        parser.add_argument('--plot', dest='plot', action='store_true', default=False)
        parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='output location, default: std out', default=sys.stdout)

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return ReadCountAnalysis(simArgs)

class Evidences(object):
    pass

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

        self.features = HTSeq.GenomicArrayOfSets("auto", stranded=True)

        for feature in self.genomeAnnotation:
            if feature.type == "gene":
                self.features[feature.iv] += self.getFeatureID(feature)

    def readReadInfo(self, args):

        self.readInfo = {}

        if args.read_info != None:
            for readInfoFile in args.read_info:
                allReadData = DataFrame.parseFromFile(readInfoFile, cDelim='\t')

                for row in allReadData:
                    readName = row["READ_NAME"]
                    readType = row["TYPE"]
                    readRun = row["USER_RUN_NAME"]
                    readLength = int(row["READ_LENGTH"])

                    self.readInfo[readName] = (readType, readRun, readLength)

    def prepareInputs(self, args):
        self.readReadInfo(args)
        self.readGenomeAnnotation(args)

        # TODO fix write out not overwriting header ...
        self.writeLinesToOutput(args.output, "\t".join(['gene', 'coverage', 'coverage_rank', 'primary_read_count', 'primary_read_rank', 'read_counts', 'read_counts_rank']) + "\n", mode='w')


        return args.sam

    def getFeatureID(self, feature):
        featureLocusName = feature.attr['gene_id'] if 'gene_id' in feature.attr else feature.name
        return featureLocusName

    def execParallel(self, data, environment):

        cvg = HTSeq.GenomicArray("auto", stranded=False, typecode='i')

        foundReadsAligned = set()
        foundReadsNotAligned = set()

        readCounts = Counter()
        readCountsPrimary = Counter()
        alignmentCount = 0
        alignmentsPerRead = Counter()

        for file in data:
            alignment_file = HTSeq.SAM_Reader( file )



            for alngt in alignment_file:

                if not alngt.read.name in self.readInfo:
                    self.readInfo[alngt.read.name] = (Fast5TYPE.UNKNOWN, None, len(alngt.read.seq))


                if alngt.aligned:

                    alignmentCount += 1

                    alignmentsPerRead[alngt.read.name] += 1

                    foundReadsAligned.add( alngt.read.name )
                    cigars = alngt.cigar
                    for cigar in cigars:

                        if cigar.type == 'M':
                            cvg[cigar.ref_iv] += 1

                    gene_ids = set()
                    for iv, val in self.features[alngt.iv].steps():
                        gene_ids |= val
                    if len(gene_ids) == 1:
                        gene_id = list(gene_ids)[0]
                        readCounts[gene_id] += 1

                        if alngt.not_primary_alignment == False:
                            readCountsPrimary[gene_id] += 1

                    elif len(gene_ids) == 0:
                        pass
                    else:

                        for gene_id in gene_ids:
                            readCounts[gene_id] += 1

                        if alngt.not_primary_alignment == False:
                            for gene_id in gene_ids:
                                readCountsPrimary[gene_id] += 1


                else:
                    foundReadsNotAligned.add( alngt.read.name )


        readsPerAlignment = Counter()
        for x in alignmentsPerRead:
            ac = alignmentsPerRead[x]
            readsPerAlignment[ac] += 1

        for ac in sorted([x for x in readsPerAlignment]):
            print("Reads with", ac, "alignments: ", readsPerAlignment[ac])

        print("Total alignments", alignmentCount)
        print("Found aligned reads: " + str(len(foundReadsAligned)))
        print("Found unaligned reads: " + str(len(foundReadsNotAligned)))
        print("Total reads: " + str(len(foundReadsAligned) + len(foundReadsNotAligned)))

        totalBaseCount = sum([ self.readInfo[x][2] for x in self.readInfo  ])
        print("Total bases in files: " + str(totalBaseCount))

        featureLengths = {}
        featureCoverage = {}
        for feature in self.genomeAnnotation:

            if feature.type == "gene":

                featureLocusName = self.getFeatureID(feature)

                if not featureLocusName in featureLengths:
                    featureLengths[ featureLocusName ] = []

                featureLengths[ featureLocusName ].append(feature.iv.end - feature.iv.start)

                if not feature.name in featureCoverage:
                    featureCoverage[ featureLocusName ] = []

                featureCoverage[ featureLocusName ].append(cvg[feature.iv])

        allKeys = []

        for x in featureCoverage:
            allKeys.append(x)


        coverages = defaultdict(int)

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
                coverages[x] = cov


        allGeneNames = set()
        for x in coverages:
            allGeneNames.add(x)
        for x in readCounts:
            allGeneNames.add(x)
        for x in readCountsPrimary:
            allGeneNames.add(x)

        allGeneNames = list(allGeneNames)
        rankedCoverage = stats.rankdata( [ coverages[gene] for gene in allGeneNames ] )
        rankedReadCount = stats.rankdata( [ readCounts[gene] for gene in allGeneNames] )
        rankedReadCountsPrimary = stats.rankdata([readCountsPrimary[gene] for gene in allGeneNames])

        allRankedCovs = []
        for i in range(0, len(allGeneNames)):

            elem = Evidences()
            elem.name = allGeneNames[i]

            elem.coverage = coverages[elem.name]
            elem.coverage_rank = rankedCoverage[i]

            elem.read_count = readCounts[elem.name]
            elem.read_rank = rankedReadCount[i]

            elem.primary_read_rank = rankedReadCountsPrimary[i]
            elem.primary_read_count = readCountsPrimary[elem.name]

            allRankedCovs.append( elem )


        allLines = []
        for x in sorted(allRankedCovs, key=lambda x: x.read_count):

            dataStr = "\t".join([ x.name, str(x.coverage), str(x.coverage_rank), str(x.primary_read_count), str(x.primary_read_rank), str(x.read_count), str(x.read_rank) ])
            allLines.append( dataStr + "\n" )

        self.writeLinesToOutput(environment.output, allLines)

        return allRankedCovs


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = []

        return existResult + newResult


    def makeResults(self, parallelResult, oEnvironment, args):

        if args.plot:

            plotData = {}
            plotData['coverage'] = [x.coverage for x in parallelResult]
            PorePlot.plotCumHistogram( plotData, [x for x in plotData], 'Coverage Plot', bins=len(parallelResult), pltcfg=args.pltcfg )

            plotData = {}
            plotData['read_count'] = [x.read_count for x in parallelResult]
            PorePlot.plotCumHistogram(plotData, [x for x in plotData], 'Read Count Plot', bins=len(parallelResult),
                                      pltcfg=args.pltcfg)

            plotData = Counter()
            for x in parallelResult:
                plotData[x.name] = x.read_count

            doPlotData = {"Read Count": {}}
            for x in plotData.most_common(20):
                doPlotData["Read Count"][x[0]] =  x[1]
            PorePlot.plotBars(doPlotData, "Read Counts", "xlabel", "ylabel", xlabelrotation="vertical", pltcfg=args.pltcfg)






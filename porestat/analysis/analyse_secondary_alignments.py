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


class SecondaryAlignmentAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(SecondaryAlignmentAnalysisFactory, self).__init__(parser, self._addParser(subparsers, which), which)

    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-s', '--sam', nargs='+', type=argparse.FileType('r'), required=True, help='alignment files')
        parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='output location, default: std out', default=sys.stdout)
        parser.add_argument('-a', '--all', default=False, action="store_true")
        parser.add_argument('-f', '--fasta', type=argparse.FileType('r'), help='de table file to read in',
                            required=True)

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return SecondaryAlignmentAnalysis(simArgs)

class Evidences(object):
    pass

class SecondaryAlignmentAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(SecondaryAlignmentAnalysis, self).__init__(args)

        self.readInfo = None
        self.genomeAnnotation = None

    def _makePropDict(self):

        props = {}

        props['COUNTS'] = {}
        props['COUNTS']['ALL'] = 0

        return props

    def readGenomeAnnotation(self, args):

        pass

    def readReadInfo(self, args):

        pass

    def prepareInputs(self, args):
        self.readReadInfo(args)
        self.readGenomeAnnotation(args)

        # TODO fix write out not overwriting header ...
        self.writeLinesToOutput(args.output, "\t".join(['read.id', 'target_chrom', 'target_start', 'target_end', 'reverse_align', 'as']) + "\n", mode='w')

        return args.sam

    def getFeatureID(self, feature):
        featureLocusName = feature.attr['gene_id'] if 'gene_id' in feature.attr else feature.name
        return featureLocusName

    def gilength(self, iv1):
        return iv1.end-iv1.start+1

    def findOptionalField(self, alngt, field, default=None):

        for key,val in alngt.optional_fields:

            if key.upper() == field.upper():
                return val

        return default

    def execParallel(self, data, environment):

        allSeqs = dict( (s[1], str(s[0]).upper()) for s in HTSeq.FastaReader(self.args.fasta.name, raw_iterator=True) )

        algntNames = defaultdict(list)

        gapOpenC = -5
        gapExtendC = -2
        matchC = 1
        mismatchC = -4

        def scoreGap( gaplength ):

            if gaplength == 1:
                return mismatchC

            return gapOpenC + gapExtendC * gaplength

        def scoreMatch( seq1, seq2):

            score = 0
            for ci, cj in zip(seq1, seq2):
                if ci == cj:
                    score += matchC
                else:
                    score += mismatchC

            return score

        read2seq = {}

        for file in data:
            alignment_file = HTSeq.SAM_Reader( file )

            for alngt in alignment_file:

                alngtName = None

                if alngt.read != None:
                    alngtName = alngt.read.name

                if alngtName != None:
                    alngtName = alngtName.split(" ")[0]




                if alngt.aligned:

                    alngtChrom = None

                    if alngt.iv != None:
                        alngtChrom = alngt.iv.chrom

                    if alngtChrom != None:
                        alngtChrom = alngtChrom.split(" ")[0]

                    if alngt.read.seq != b'*':
                        read2seq[alngtName] = alngt.read

                    #if alngt.supplementary:
                    #    continue

                    if alngt.read.seq == b'*' or alngt.read.seq == b'':
                        readstring = read2seq.get(alngtName, "")
                    else:
                        readstring = alngt.read

                    if alngt.flag & 0x10:
                        readstring = readstring.get_reverse_complement()

                    readstring = readstring.seq.decode()

                    lastGapIdx = len(alngt.cigar)
                    for idx, cigarOp in enumerate(reversed(alngt.cigar)):
                        if cigarOp.type in ['I', 'S', 'H', 'D']:
                            lastGapIdx = len(alngt.cigar)-1-idx
                            break

                    readScoreGlobal = 0
                    readScoreFreeshift = 0
                    firstGap = True

                    for idx, cigarOp in enumerate(alngt.cigar):

                        if cigarOp.type in ['I', 'S', 'H']:
                            readScoreGlobal += scoreGap( cigarOp.query_to-cigarOp.query_from )
                            if not firstGap and idx < lastGapIdx:
                                readScoreFreeshift += scoreGap( cigarOp.ref_iv.end-cigarOp.ref_iv.start )

                            firstGap = False
                        elif cigarOp.type in ['D']:
                            readScoreGlobal += scoreGap( cigarOp.ref_iv.end-cigarOp.ref_iv.start )

                            if not firstGap and idx < lastGapIdx:
                                readScoreFreeshift += scoreGap( cigarOp.ref_iv.end-cigarOp.ref_iv.start )

                            firstGap = False

                        elif cigarOp.type in ["M"]:

                            readCigar = readstring[cigarOp.query_from:cigarOp.query_to]
                            refCigar = allSeqs[alngtChrom][cigarOp.ref_iv.start:cigarOp.ref_iv.end]

                            matchScore = scoreMatch( readCigar, refCigar)
                            readScoreGlobal += matchScore
                            readScoreFreeshift += matchScore
                            firstGap = False

                        elif cigarOp.type in ["N"]:
                            continue
                        else:
                            print("Unknown Cigar", cigarOp)

                    readScoreGlobalByLength = readScoreGlobal / (alngt.iv.end-alngt.iv.start+1)
                    readScoreFreeshiftByLength = readScoreFreeshift / (alngt.iv.end - alngt.iv.start + 1)

                    adata = (
                        alngtChrom,
                        alngt.iv.start,
                        alngt.iv.end,
                        alngt.iv.strand,
                        self.findOptionalField(alngt, "AS", -1),
                        readScoreGlobal,
                        readScoreFreeshift,
                        readScoreGlobalByLength,
                        readScoreFreeshiftByLength,
                        readstring)
                    algntNames[alngtName].append(
                        adata
                    )

        print("Writing out files")

        lines2write = []
        for readName in algntNames:

            if len(algntNames[readName]) > 1 or self.args.all:

                for readdata in algntNames[readName]:
                    dataStr = "\t".join(
                        [str(x) for x in [readName] + list(readdata)]
                    ) + "\n"

                    lines2write.append(dataStr)

        self.writeLinesToOutput(environment.output, lines2write)


        return []


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = []

        return existResult + newResult


    def makeResults(self, parallelResult, oEnvironment, args):

        pass


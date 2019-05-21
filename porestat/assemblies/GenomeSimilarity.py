import argparse
import itertools
import pysam

import sys
import tempfile

import numpy
from Bio.Blast import NCBIXML
from matplotlib import gridspec
from porestat.hdf5tool import FASTQ
from porestat.tools.PTToolInterface import PSToolInterfaceFactory, PSToolInterface
from Bio import SeqIO
from collections import OrderedDict, defaultdict

from matplotlib import pyplot as plt

import mappy as mp


class GenomeSimilarityPlotFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(GenomeSimilarityPlotFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-fq', '--fastq', type=argparse.FileType('r'), help='number of reads to extract', required=True)
        parser.add_argument('-fq2', '--fastq2', type=argparse.FileType('r'), help='number of reads to extract',
                            required=False, default=None)
        parser.add_argument('-a', '--seqA', type=argparse.FileType('r'), help='name of file with sequences (reference)', required=True)
        parser.add_argument('-b', '--seqB', type=argparse.FileType('r'), help='name of file with sequences (query)', required=True)
        parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='output image', required=True)
        parser.add_argument('-t', '--tmp', type=str, help='path to tmp folder', default=tempfile.gettempdir())
        parser.add_argument('-ill', '--short', action="store_true", default=False, help="short illumina reads")

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return GenomeSimilarityPlotter(simArgs)

import os

class Environment(object):
    pass


class GenomeSimilarityPlotter(PSToolInterface):

    def __init__(self, args):

        super(GenomeSimilarityPlotter, self).__init__(args)

    def getFileName(self, fname):
        return os.path.join(self.args.tmp, fname)

    def prepareDatabase(self, inFile):

        dbOutFile = self.getFileName("bdb")
        dbcall = "makeblastdb -in {infile} -dbtype nucl -input_type fasta -out {ofile}".format(infile=inFile.name, ofile=dbOutFile)

        os.system(dbcall)

        return dbOutFile

    def runAgainstDatabase(self, dbFilePref, querySeq):

        oname = self.getFileName("query_res.xml")
        blastcall = "blastn -query {qfname} -db {db} -outfmt 5 -out {outname}".format(qfname=querySeq.name, db=dbFilePref, outname=oname)
        os.system(blastcall)

        print(blastcall)

        return oname

    def makeQuery(self):
        pass

    def prepareDotPlot(self):

        dbFile = self.prepareDatabase(self.args.seqA)

        alignmentFile = self.runAgainstDatabase(dbFile, self.args.seqB)

        blast_record = None

        with open( alignmentFile, 'r') as fin:
            blast_record = list(NCBIXML.parse(fin))

        return blast_record

    def acceptReadsBam(self, hits, hasRead1, hasRead2, hasReads2):

        if len(hits) == 0:
            return False

        exHit = hits[0]

        isPaired = exHit.is_paired

        hitsperread = {1: 0, 2: 0}

        for hit in hits:
            hitsperread[1 if hit.is_read1 else 2] += 1

        if not isPaired:
            if hitsperread[1] == 1:
                return True

        else:

            if isPaired:
                # read1 and read2 must have 1 alignment
                if hitsperread[1] == 1 and hitsperread[2] == 1:
                    return True

            elif 1 in hitsperread and not 2 in hitsperread:

                if hitsperread[1] == 1:
                    return True

            elif not 1 in hitsperread and 2 in hitsperread:

                if hitsperread[2] == 1:
                    return True

        return False


    def getMate(self, read, hits):

        for x in hits:

            if (x.is_read2 and read.is_read1) or (x.is_read1 and read.is_read2):

                if read.next_reference_id == x.reference_id and read.next_reference_start == x.reference_start:

                    return x

        return None



    def prepareCoverageFromBam(self, bamFile):


        bamFile = pysam.AlignmentFile(bamFile)

        refs = bamFile.nreferences

        seqNames = []
        covArrayPrim = {}
        covArraySec = {}
        for i in range(0, refs):

            seqName = bamFile.get_reference_name(i)
            seqLenght = bamFile.get_reference_length(i)

            seqNames.append(seqName)

            covArrayPrim[seqName] = numpy.zeros(seqLenght, 1)
            covArraySec[seqName] = numpy.zeros(seqLenght, 1)


        curHits = []
        curReadName = None

        dupHits = {}

        for algn in bamFile:

            if curReadName != algn.query_name:

                # check alignments
                if self.acceptReadsBam(curHits):

                    for read in curHits:
                        if read.is_read1:

                            read2 = self.getMate(read, curHits)

                            if read2 == None:
                                continue

                            for q,r in read.get_aligned_pairs(matches_only=True):
                                covArrayPrim[read.reference_name][r, 0] += 1

                            for q,r in read2.get_aligned_pairs(matches_only=True):
                                covArrayPrim[read2.reference_name][r, 0] += 1

                else:

                    if len(curHits) > 0:
                        dupHits[curReadName] = curHits

                curHits = []
                curReadName = algn.query_name

            curHits.append(algn)

        for rname in dupHits:

            curHits = dupHits[rname]










    def prepareCoverage(self, seqFile):

        preset = "map-ont"

        if self.args.short or self.args.fastq2 != None:
            preset = "sr"
            print("Short Reads Mode")
            print(self.args.fastq.name)
            print(self.args.fastq2.name)

        a = mp.Aligner(seqFile.name, n_threads=2, preset=preset)  # load or build index
        if not a:
            raise Exception("ERROR: failed to load/build index")


        seqNames = a.seq_names
        covArrayPrim = {}
        covArraySec = {}

        for name in seqNames:
            seq = a.seq(name, start=0, end=0x7fffffff)
            covArrayPrim[name] = numpy.zeros((len(seq), 1))
            covArraySec[name] = numpy.zeros((len(seq), 1))


        maxAlignments = -1

        alignedReads = 0

        print("Alignment process 1")
        aligned = 0

        reads1 = mp.fastx_read(self.args.fastq.name)
        reads2 = [None]
        hasReads2 = False

        if self.args.fastq2 != None:
            reads2 = mp.fastx_read(self.args.fastq2.name)
            hasReads2 = True

        readNamesStage2 = set()

        for read1, read2 in itertools.zip_longest(reads1, reads2):

            if maxAlignments >= 0 and aligned > maxAlignments:
                break

            if read1 == None and read2 == None:
                continue

            if read1 != None:
                name1, seq1, qual1 = read1
            else:
                seq1 = None
                name1 = None

            if read2 != None:
                name2, seq2, qual2 = read2
            else:
                seq2 = None
                name2 = None

            hits = [x for x in a.map(seq1, seq2, cs=False, MD=False)]

            if not self.acceptReads(hits, read1 != None, read2 != None, hasReads2):
                readNamesStage2.add( (name1, name2) )
                continue

            alignedReads += 1

            for hit in hits:
                if not hit.is_primary:
                    continue

                refSeqName = hit.ctg

                aligned += 1

                for i in range(hit.r_st-1, hit.r_en):
                    covArrayPrim[refSeqName][i,0] += 1


        covArrayMerged = {}
        for x in covArrayPrim:
            covArrayMerged[x] = numpy.copy(covArrayPrim[x])

        print("Alignment process 2")
        aligned = 0

        reads1 = mp.fastx_read(self.args.fastq.name)
        reads2 = [None]

        if self.args.fastq2 != None:
            reads2 = mp.fastx_read(self.args.fastq2.name)

        for read1, read2 in itertools.zip_longest(reads1, reads2):

            if maxAlignments >= 0 and aligned > maxAlignments:
                break

            if read1 == None and read2 == None:
                continue

            if read1 != None:
                name1, seq1, qual1 = read1
            else:
                seq1 = None
                name1 = None

            if read2 != None:
                name2, seq2, qual2 = read2
            else:
                seq2 = None
                name2 = None

            if not (name1, name2) in readNamesStage2:
                continue

            hits = [x for x in a.map(seq1, seq2, cs=False, MD=False)]

            position2cov = []

            for hit in hits:

                refSeqName = hit.ctg

                hitPosition = (refSeqName, hit.r_st, hit.r_en, self.getCoverageForRegion(covArrayMerged[refSeqName], hit))
                position2cov.append(hitPosition)

                if hit.is_primary:
                    for i in range(hit.r_st-1, hit.r_en):
                        covArraySec[refSeqName][i,0] += 1

            if len(position2cov) == 0:
                continue

            alignedReads += 1


            hitPositions = sorted(position2cov, key=lambda x: x[3])
            hitPosition = hitPositions[0]
            for i in range(hitPosition[1]-1, hitPosition[2]):
                covArrayMerged[hitPosition[0]][i,0] += 1

        return (covArrayPrim, covArraySec, covArrayMerged, alignedReads)

    def getCoverageForRegion(self, covArrayMerged, hit):

        subArray = covArrayMerged[(hit.r_st-1):hit.r_en, 0]
        return numpy.mean(subArray)

    def acceptReads(self, hits, hasRead1, hasRead2, hasReads2):

        hitsperread = {1: 0, 2: 0}

        for hit in hits:
            hitsperread[hit.read_num] += 1

        if not hasReads2:
            if hitsperread[1] == 1:
                return True

        else:

            if hasRead1 and hasRead2:
                # read1 and read2 must have 1 alignment
                if hitsperread[1] == 1 and hitsperread[2] == 1:
                    return True

            elif hasRead1 and not hasRead2:

                if hitsperread[1] == 1:
                    return True

            elif not hasRead1 and hasRead2:

                if hitsperread[2] == 1:
                    return True


        return False


    def plotAlignment(self, ax, seqA, seqB, blast):

        alignRegions = []

        for blastRes in blast:

            querySeq = blastRes.query

            if not querySeq == seqA:
                continue

            for alignment in blastRes.alignments:

                subjSeq = alignment.hit_def

                if not subjSeq == seqB:
                    continue

                for hsp in alignment.hsps:

                    if hsp.align_length < 1000:
                        continue

                    if hsp.expect > 0.05:
                        continue

                    alignPart = (
                    hsp.query_start,
                    hsp.query_end,
                    hsp.sbjct_start,
                    hsp.sbjct_end
                    )

                    alignRegions.append(alignPart)



        for alignPart in alignRegions:

            ax.plot([alignPart[0], alignPart[1]],[alignPart[2], alignPart[3]])

        ax.set_xlabel(seqA)
        ax.set_ylabel(seqB)

        return ax







    def exec(self):

        aSeqs = {}
        bSeqs = {}

        for name, seq, _ in mp.fastx_read(self.args.seqA.name):
            aSeqs[name] = seq

        for name, seq, _ in mp.fastx_read(self.args.seqB.name):
            bSeqs[name] = seq

        print([x for x in aSeqs])
        print([x for x in bSeqs])


        covDataPrimA, covDataSecA, covDataMergedA, alignedA = self.prepareCoverage(self.args.seqA)
        covDataPrimB, covDataSecB, covDataMergedB, alignedB = self.prepareCoverage(self.args.seqB)

        alignments = self.prepareDotPlot()

        fig = plt.figure(figsize=(8,8))

        gs = gridspec.GridSpec(len(aSeqs)+1, len(bSeqs)+1, width_ratios=[1] + [3]*len(bSeqs), height_ratios=[3]*len(aSeqs) + [1] )

        for aIdx, aSeq in enumerate(aSeqs):

            ax = fig.add_subplot(gs[aIdx*(len(bSeqs)+1)])

            aY = covDataPrimA[aSeq]
            ax.plot(aY, numpy.arange(0, aY.shape[0]))

            aY = covDataSecA[aSeq]
            ax.plot(aY,numpy.arange(0, aY.shape[0]))

            aY = covDataMergedA[aSeq]
            ax.plot(aY,numpy.arange(0, aY.shape[0]))

            ax.invert_xaxis()
            ax.set_xlabel("Aligned Reads = {nr}".format(nr=alignedA))

            for bIdx, bSeq in enumerate(bSeqs):
                ax = fig.add_subplot(gs[aIdx * (len(bSeqs) + 1) + 1 + bIdx])

                self.plotAlignment(ax, bSeq, aSeq, alignments)


        for bIdx, bSeq in enumerate(bSeqs):
            ax = fig.add_subplot(gs[len(aSeqs) * (len(bSeqs) + 1) + 1 + bIdx])
            bY = covDataPrimB[bSeq]
            ax.plot(numpy.arange(0, bY.shape[0]), bY)

            bY = covDataSecB[bSeq]
            ax.plot(numpy.arange(0, bY.shape[0]), bY)

            bY = covDataMergedB[bSeq]
            ax.plot(numpy.arange(0, bY.shape[0]), bY)

            ax.set_xlabel("Aligned Reads = {nr}".format(nr=alignedB))

        plt.savefig(self.args.output.name + ".pdf", bbox_inches="tight")
        plt.savefig(self.args.output.name + ".png", bbox_inches="tight")

        plt.tight_layout()
        plt.show()

        return True

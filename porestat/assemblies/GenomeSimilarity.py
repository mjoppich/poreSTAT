import argparse
import datetime
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
        parser.add_argument('-fq', '--fastq', type=argparse.FileType('r'), help='number of reads to extract', required=False, default=None)
        parser.add_argument('-fq2', '--fastq2', type=argparse.FileType('r'), help='number of reads to extract',
                            required=False, default=None)

        parser.add_argument('-ba', '--bamA', type=argparse.FileType('rb'), help='name of file with sequences (reference)',
                            required=False, default=None)
        parser.add_argument('-bb', '--bamB', type=argparse.FileType('rb'), help='name of file with sequences (query)',
                            required=False, default=None)

        parser.add_argument('-a', '--seqA', type=argparse.FileType('r'), help='name of file with sequences (reference)', required=True)
        parser.add_argument('-b', '--seqB', type=argparse.FileType('r'), help='name of file with sequences (query)', required=True)
        parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='output image', required=True)
        parser.add_argument('-t', '--tmp', type=str, help='path to tmp folder', default=tempfile.gettempdir())
        parser.add_argument('-ill', '--short', action="store_true", default=False, help="short illumina reads")

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        parser= argparse.ArgumentParser()

        if args.fastq2 != None and not args.fastq != None:
            parser.error("--fastq2 given, but not --fastq")

        if args.fastq != None and (args.bamA != None or args.bamB != None):
            parser.error('if --fastq given, you may not give --bamA or --bamB')

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

    def acceptReadsBam(self, hits):

        hitsLen = len(hits)
        if hitsLen == 0:
            return False

        exHit = hits[0]
        isPaired = exHit.is_paired

        if isPaired and hitsLen == 1:
            return False

        if not isPaired and hitsLen == 1:
            return True

        if hitsLen == 2:
            if isPaired:
                if hits[0].is_read1 and hits[1].is_read2:
                    return True
                elif hits[1].is_read1 and hits[0].is_read2:
                    return True
                else:
                    return False

        if hitsLen > 2 and isPaired:
            return False


        hitsperread = {1: 0, 2: 0}

        for hit in hits:
            hitsperread[1 if hit.is_read1 else 2] += 1

        if not isPaired:
            if hitsperread[1] == 1:
                return True

        elif isPaired:
            # read1 and read2 must have 1 alignment
            if hitsperread[1] == 1 and hitsperread[2] == 1:
                return True
            else:
                return False

        return False


    def getMate(self, read, hits):

        readNextRefID = read.next_reference_id
        readNextRefStart = read.next_reference_start
        readRead1 = read.is_read1

        for x in hits:

            if x.is_read1 == readRead1:
                continue

            if readNextRefID == x.reference_id and readNextRefStart == x.reference_start:
                return True, x

        return False, None

    def bamHandleHits(self, curHits, covArray):

        alignedReads = 0
        for read in curHits:
            if read.is_read1:

                ret, read2 = self.getMate(read, curHits)

                if ret == False:
                    continue

                alignedReads += 2

                refName = read.reference_name
                refCovArray = covArray[refName]

                for r in range(read.reference_start, read.reference_end):
                    refCovArray[0,r] += 1

                refName = read2.reference_name
                refCovArray = covArray[refName]

                for r in range(read2.reference_start, read2.reference_end):
                    refCovArray[0,r] += 1

        return alignedReads

    def processAlignmentsBam(self, dupHits, covArray):

        alignedReads = 0
        namestoDelete = []
        print(datetime.datetime.now(), "Checking", len(dupHits), "Read Names")

        for rname in dupHits:
            thits = dupHits[rname]

            if self.acceptReadsBam(thits):
                alignedR = self.bamHandleHits(thits, covArray)
                alignedReads += alignedR

                if alignedR > 0:
                    namestoDelete.append(rname)

        print(datetime.datetime.now(), "Merged", len(namestoDelete), "Read Names")
        for rname in namestoDelete:
            del dupHits[rname]

        return alignedReads

    def prepareCoverageFromBam(self, bamFile):


        bamFile = pysam.AlignmentFile(bamFile.name)

        refs = bamFile.nreferences

        seqNames = []
        covArrayPrim = {}
        covArraySec = {}
        covArrayClipped = {}

        for i in range(0, refs):

            seqName = bamFile.get_reference_name(i)
            seqLenght = bamFile.get_reference_length(seqName)

            seqNames.append(seqName)

            covArrayPrim[seqName] = numpy.zeros((1,seqLenght))
            covArraySec[seqName] = numpy.zeros((1,seqLenght))
            covArrayClipped[seqName] = numpy.zeros((1,seqLenght))


        curHits = []
        curReadName = None
        alignedReads = 0
        dupHits = {}

        maxReads = -1

        print(datetime.datetime.now(), "Processing Bam File")
        print(datetime.datetime.now(), "Accept 1")

        for aidx, algn in enumerate(bamFile):

            if aidx % 100000 == 0:
                print(datetime.datetime.now(), "Processing alignment", aidx)

                alignedReads += self.processAlignmentsBam(dupHits, covArrayPrim)



            if algn.is_unmapped:
                continue

            if maxReads != -1 and alignedReads > maxReads:
                break

            if curReadName != algn.query_name:
                if len(curHits) > 0:
                    if curReadName in dupHits:
                        dupHits[curReadName] += curHits
                    else:
                        dupHits[curReadName] = curHits

                curHits = []
                curReadName = algn.query_name

            curHits.append(algn)

        # process whatever is left over
        alignedReads += self.processAlignmentsBam(dupHits, covArrayPrim)

        covArrayMerged = {}
        for x in covArrayPrim:
            covArrayMerged[x] = numpy.copy(covArrayPrim[x])

        print(datetime.datetime.now(), "Accept 2")
        print(datetime.datetime.now(), "Reprocessing elements:", len(dupHits))


        zeroPos2Cov = 0
        alignedReads = 0
        for rname in dupHits:

            if maxReads != -1 and alignedReads > maxReads:
                break

            curHits = dupHits[rname]
            position2cov = []

            for read in curHits:

                if read.is_read1:

                    ret, read2 = self.getMate(read, curHits)
                    refSeqName = read.reference_name

                    cov1 = self.getMeanCoverage(covArrayMerged[refSeqName], read.reference_start+1, read.reference_end+1)

                    if ret:
                        cov2 = self.getMeanCoverage(covArrayMerged[refSeqName], read2.reference_start+1, read2.reference_end+1)
                    else:
                        cov2 = 0

                    rcount = 1

                    if ret == True:
                        rcount += 1

                    hitPosition = (
                        refSeqName, read.reference_start, (cov1+cov2)/rcount, rcount, read, read2
                    )
                    position2cov.append(hitPosition)

                    alignedReads += rcount

                    if not read.is_secondary and not read.is_supplementary:
                        for i in range(read.reference_start, read.reference_end):
                            covArraySec[refSeqName][0,i] += 1

            if len(position2cov) == 0:
                zeroPos2Cov += 1
                continue

            hitPositions = sorted(position2cov, key=lambda x: x[2])
            hitPosition = hitPositions[0]

            refName = hitPosition[4].reference_name
            refCovArray = covArrayMerged[refName]

            for r in range(hitPosition[4].reference_start, hitPosition[4].reference_end):
                refCovArray[0, r] += 1

            if hitPosition[3] == 2:
                refName = hitPosition[5].reference_name
                refCovArray = covArrayMerged[refName]

                for r in range(hitPosition[5].reference_start, hitPosition[5].reference_end):
                    refCovArray[0, r] += 1

        print(datetime.datetime.now(), "Zero Pos2 Cov", zeroPos2Cov)

        return (covArrayPrim, covArraySec, covArrayMerged, covArrayClipped, alignedReads)


    def addHitCoverage(self, hit, refSeqName, covArray):

        for i in range(hit.r_st - 1, hit.r_en):
            covArray[refSeqName][0,i] += 1

    def handleClippedCigar(self, aligner, cseq, cigarElem, refSeqName, covArray):

        #softClipped = 4
        #hardClipped = 5

        if cigarElem[1] in [4,5]:
            if cigarElem[0] > 500:
                nseq = cseq[0:cigarElem[0]]
                nhits = [x for x in aligner.map(nseq, cs=False, MD=False)]

                for nhit in nhits:
                    if nhit.is_primary:
                        self.addHitCoverage(nhit, refSeqName, covArray)

                return 1
        return 0

    def prepareCoverage(self, seqFile):

        preset = None#"map-ont"

        if self.args.short or self.args.fastq2 != None:
            preset = "sr"
            print("Short Reads Mode")
            print(self.args.fastq.name)
            print(self.args.fastq2.name)

        a = mp.Aligner(seqFile.name, n_threads=2, preset=preset, )  # load or build index
        if not a:
            raise Exception("ERROR: failed to load/build index")


        seqNames = a.seq_names
        covArrayPrim = {}
        covArraySec = {}
        covArrayClipped = {}

        for name in seqNames:
            seq = a.seq(name, start=0, end=0x7fffffff)
            covArrayPrim[name] = numpy.zeros((1,len(seq)))
            covArraySec[name] = numpy.zeros((1,len(seq)))
            covArrayClipped[name] = numpy.zeros((1,len(seq)))

        maxAlignments = -1
        alignedReads = 0
        ccount = 0

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

                self.addHitCoverage(hit, refSeqName, covArrayPrim)
                # check for clipped parts

                if hit.read_num == 1:
                    cseq = seq1
                else:
                    cseq = seq2

                startCigar = hit.cigar[0]
                ccount += self.handleClippedCigar(a, cseq, startCigar, refSeqName, covArrayClipped)
                endCigar = hit.cigar[-1]
                ccount += self.handleClippedCigar(a, cseq, endCigar, refSeqName, covArrayClipped)




        covArrayMerged = {}
        for x in covArrayPrim:
            covArrayMerged[x] = numpy.copy(covArrayPrim[x])

        print("Clipped Proc 1", ccount)
        print("Alignment process 2")
        aligned = 0

        reads1 = mp.fastx_read(self.args.fastq.name)
        reads2 = [None]

        if self.args.fastq2 != None:
            reads2 = mp.fastx_read(self.args.fastq2.name)

        zeroPos2Cov = 0
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

                hitPosition = (refSeqName,
                               hit.r_st,
                               self.getCoverageForRegion(covArrayMerged[refSeqName], hit),
                               hit
                               )
                position2cov.append(hitPosition)

                if hit.is_primary:
                    for i in range(hit.r_st-1, hit.r_en):
                        covArraySec[refSeqName][0,i] += 1

            if len(position2cov) == 0:
                zeroPos2Cov += 1

                continue

            alignedReads += 1


            hitPositions = sorted(position2cov, key=lambda x: x[2])
            hitPosition = hitPositions[0]
            self.addHitCoverage(hitPosition[3], hitPosition[0], covArrayMerged)
            # check for clipped parts

            if hitPosition[3].read_num == 1:
                cseq = seq1
            else:
                cseq = seq2

            startCigar = hitPosition[3].cigar[0]
            ccount += self.handleClippedCigar(a, cseq, startCigar, hitPosition[0], covArrayClipped)
            endCigar = hitPosition[3].cigar[-1]
            ccount += self.handleClippedCigar(a, cseq, endCigar, hitPosition[0], covArrayClipped)

        print("0 cov 2 cov", zeroPos2Cov)
        print("Clipped Proc 2", ccount)
        return (covArrayPrim, covArraySec, covArrayMerged, covArrayClipped, alignedReads)

    def getCoverageForRegion(self, covArray, hit):

        if hit == None:
            return 0

        return self.getMeanCoverage(covArray, hit.r_st, hit.r_en)

    def getMeanCoverage(self, covArray, start, end):

        subArray = covArray[(start-1):end, 0]
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




    def plotCovData(self, ax, seqName, covDataPrim, covDataSec, covDataMerged, covDataClipped, alignedCount, vertical=False):

        pl = None
        sl = None
        ml = None
        cl = None

        if covDataPrim != None:
            aY = covDataPrim[seqName][0]
            rangeElems = numpy.arange(0, len(aY))

            if vertical:
                pl, = ax.plot(aY, rangeElems, color="blue")
            else:
                pl, = ax.plot(rangeElems, aY, color="blue")

        if covDataSec != None:
            aY = covDataSec[seqName][0]
            rangeElems = numpy.arange(0, len(aY))

            if vertical:
                sl, = ax.plot(aY, rangeElems, color="orange")
            else:
                sl, = ax.plot(rangeElems, aY, color="orange")

        if covDataMerged != None:
            aY = covDataMerged[seqName][0]
            rangeElems = numpy.arange(0, len(aY))

            if vertical:
                ml, = ax.plot(aY, rangeElems, color="green")
            else:
                ml, = ax.plot(rangeElems, aY, color="green")

        if covDataClipped != None:
            aY = covDataClipped[seqName][0]
            rangeElems = numpy.arange(0, len(aY))

            if vertical:
                cl, = ax.plot(aY, rangeElems, color="red")
            else:
                cl, = ax.plot(rangeElems, aY, color="red")

        if vertical:
            ax.invert_xaxis()

        ax.set_xlabel("Aligned Reads = {nr}".format(nr=alignedCount))

        return [pl, sl, ml, cl], ["Unique", "Multiple", "Merged by Cov", "Clipped"]


    def exec(self):

        aSeqs = {}
        bSeqs = {}

        for name, seq, _ in mp.fastx_read(self.args.seqA.name):
            aSeqs[name] = seq

        for name, seq, _ in mp.fastx_read(self.args.seqB.name):
            bSeqs[name] = seq

        print([x for x in aSeqs])
        print([x for x in bSeqs])


        if self.args.fastq:
            covDataPrimA, covDataSecA, covDataMergedA, covDataClippedA, alignedA = self.prepareCoverage(self.args.seqA)
            if self.args.seqA.name == self.args.seqB.name:
                print("--seqA == --seqA, copying over results")
                covDataPrimB, covDataSecB, covDataMergedB, covDataClippedB, alignedB = covDataPrimA, covDataSecA, covDataMergedA, covDataClippedA, alignedA
            else:
                covDataPrimB, covDataSecB, covDataMergedB, covDataClippedB, alignedB = self.prepareCoverage(self.args.seqB)

        else:
            covDataPrimA, covDataSecA, covDataMergedA, covDataClippedA, alignedA = self.prepareCoverageFromBam(self.args.bamA)

            if self.args.bamA.name == self.args.bamB.name:
                print("--bamA == --bamB, copying over results")
                covDataPrimB, covDataSecB, covDataMergedB, covDataClippedB, alignedB = covDataPrimA, covDataSecA, covDataMergedA, covDataClippedA, alignedA
            else:
                covDataPrimB, covDataSecB, covDataMergedB, covDataClippedB, alignedB = self.prepareCoverageFromBam(self.args.bamB)

        alignments = self.prepareDotPlot()

        fig = plt.figure(figsize=(8,8))

        gs = gridspec.GridSpec(len(aSeqs)+1, len(bSeqs)+1, width_ratios=[1] + [3]*len(bSeqs), height_ratios=[3]*len(aSeqs) + [1] )

        for aIdx, aSeq in enumerate(aSeqs):

            ax = fig.add_subplot(gs[aIdx*(len(bSeqs)+1)])
            llines, llabels = self.plotCovData(ax, aSeq, covDataPrimA, covDataSecA, covDataMergedA, covDataClippedA, alignedA, vertical=True)

            for bIdx, bSeq in enumerate(bSeqs):
                ax = fig.add_subplot(gs[aIdx * (len(bSeqs) + 1) + 1 + bIdx])
                self.plotAlignment(ax, bSeq, aSeq, alignments)

        ax = fig.add_subplot(gs[len(aSeqs) * (len(bSeqs) + 1) ])
        ax.legend(llines, llabels)
        ax.set_axis_off()

        for bIdx, bSeq in enumerate(bSeqs):
            ax = fig.add_subplot(gs[len(aSeqs) * (len(bSeqs) + 1) + 1 + bIdx])
            self.plotCovData(ax, bSeq, covDataPrimB, covDataSecB, covDataMergedB, covDataClippedB, alignedB)


        plt.savefig(self.args.output.name + ".pdf", bbox_inches="tight")
        plt.savefig(self.args.output.name + ".png", bbox_inches="tight")

        plt.tight_layout()
        plt.show()

        return True

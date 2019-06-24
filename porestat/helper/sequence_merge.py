import argparse

import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from porestat.hdf5tool import FASTQ
from porestat.tools.PTToolInterface import PSToolInterfaceFactory, PSToolInterface
from Bio import SeqIO
from collections import OrderedDict

class FARecordSequenceMergeFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(FARecordSequenceMergeFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-f', '--fasta', nargs='+', type=argparse.FileType('r'), help='minion read folder', required=True)
        parser.add_argument('-s', '--seqs', type=int, help='minion read folder',
                            required=True)
        parser.add_argument('-n', '--name', type=str, default=None)

        def fileOpener( filename ):
            print("Opening", filename)
            open(filename, 'w').close()
            return filename

        parser.add_argument('-o', '--output', type=fileOpener, help='output location, default: std out', default=sys.stdout)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return FARecordSequenceMerge(simArgs)

import os

class Environment(object):
    pass


class FARecordSequenceMerge(PSToolInterface):

    def __init__(self, args):

        super(FARecordSequenceMerge, self).__init__(args)

    def acceptFASTQ(self, fastq):

        assert(type(fastq) == FASTQ)

        if self.args.max_length != -1:

            if len(fastq) > self.args.max_length:
                return False

        if self.args.min_length != -1:

            if len(fastq) < self.args.min_length:
                return False

        return True


    def exec(self):

        ftype = 'fasta'
        seqCount = 0

        outSeq = None
        outSeqID = self.args.name

        with open(self.args.output, 'w') as fout:
            for infile in self.args.fasta:

                seqsOfFile = 0
                for record in SeqIO.parse(infile, ftype):

                    if seqsOfFile >= self.args.seqs:
                        break

                    if outSeq == None:
                        outSeq = record.seq
                    else:
                        outSeq += record.seq

                    if outSeqID == None:
                        outSeqID = record.id


                    seqsOfFile += 1

            outRecord = SeqRecord(outSeq,
                       id=outSeqID, name=outSeqID,
                       description=outSeqID)

            print(outRecord.id)
            print(len(outRecord.seq))

            SeqIO.write([outRecord], fout, ftype)

        return True

import argparse

import sys

from porestat.hdf5tool import FASTQ
from porestat.tools.PTToolInterface import PSToolInterfaceFactory, PSToolInterface
from Bio import SeqIO
from collections import OrderedDict

class FAFQRecordExtractFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(FAFQRecordExtractFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-f', '--fasta', nargs='+', type=argparse.FileType('r'), help='minion read folder', required=True)
        parser.add_argument('-s', '--seqcount', type=int, help='number of reads to extract', required=False, default=10)
        parser.add_argument('-se', '--seqend', type=int, help='number of reads to extract', required=False, default=-1)
        parser.add_argument('-st', '--seqstart', type=int, help='number of reads to extract', required=False,
                            default=0)
        parser.add_argument('-fq', '--fastq', type=int, help='reads are fastq', default=False, required=False)
        parser.add_argument('-rv', '--reverse', action='store_true', default=False)


        def fileOpener( filename ):
            print("Opening", filename)
            open(filename, 'w').close()
            return filename

        parser.add_argument('-o', '--output', type=fileOpener, help='output location, default: std out', default=sys.stdout)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return FAFQRecordExtract(simArgs)

import os

class Environment(object):
    pass


class FAFQRecordExtract(PSToolInterface):

    def __init__(self, args):

        super(FAFQRecordExtract, self).__init__(args)

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
        if self.args.fastq:
            ftype= 'fastq'

        seqCount = 0

        with open(self.args.output, 'w') as fout:
            for infile in self.args.fasta:

                for record in SeqIO.parse(infile, ftype):

                    seqStart = 0
                    seqEnd = len(record)

                    if self.args.seqend != -1:
                        if len(record) >= self.args.seqend:
                            seqEnd = self.args.seqend

                    if self.args.seqstart != 0:
                        if self.args.seqstart > 0 and self.args.seqstart < len(record):
                            seqStart = self.args.seqstart

                    record = record[seqStart:seqEnd]

                    if self.args.reverse:
                        record = record.reverse_complement(id=True, name=True, description=True)

                    SeqIO.write([record], fout, ftype)
                    seqCount += 1

                    if seqCount >= self.args.seqcount:
                        break
                if seqCount >= self.args.seqcount:
                    break


        return True

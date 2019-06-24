import argparse

import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from porestat.hdf5tool import FASTQ
from porestat.tools.PTToolInterface import PSToolInterfaceFactory, PSToolInterface
from Bio import SeqIO
from collections import OrderedDict

class FaSeqReplaceFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(FaSeqReplaceFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-f', '--fasta', nargs='+', type=argparse.FileType('r'), help='minion read folder', required=True)

        parser.add_argument('-sstart', '--source-start', type=int, help="source start")
        parser.add_argument('-send', '--source-end', type=int, help="source end")

        parser.add_argument('-tstart', '--target-start', type=int, help="target start")
        parser.add_argument('-tend', '--target-end', type=int, help="target end")


        def fileOpener( filename ):
            print("Opening", filename)
            open(filename, 'w').close()
            return filename

        parser.add_argument('-o', '--output', type=fileOpener, help='output location, default: std out', default=sys.stdout)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return FaSeqReplace(simArgs)

import os

class Environment(object):
    pass


class FaSeqReplace(PSToolInterface):

    def __init__(self, args):

        super(FaSeqReplace, self).__init__(args)

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
        with open(self.args.output, 'w') as fout:
            revRecords = []

            for infile in self.args.fasta:

                for record in SeqIO.parse(infile, ftype):
                    print(record.id)
                    print(len(record.seq))

                    nrec = record.upper()

                    nrec.seq = record.seq[0:self.args.tstart]
                    nrec.seq += record.seq[self.args.sstart-1:self.args.send]
                    nrec.seq+= record.seq[self.args.tend:]

                    revRecords.append(nrec)



            SeqIO.write(revRecords, fout, ftype)



        return True

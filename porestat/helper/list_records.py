import argparse

import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from porestat.hdf5tool import FASTQ
from porestat.tools.PTToolInterface import PSToolInterfaceFactory, PSToolInterface
from Bio import SeqIO
from collections import OrderedDict
import re

from porestat.utils.Files import printToFile


class FaListRecordsFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(FaListRecordsFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-f', '--fasta', nargs='+', type=argparse.FileType('r'), help='minion read folder', required=True)

        def fileOpener( filename ):
            print("Opening", filename)
            open(filename, 'w').close()
            return filename

        parser.add_argument('-o', '--output', type=fileOpener, help='output location, default: std out', default=sys.stdout)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return FaListRecords(simArgs)

import os

class Environment(object):
    pass


class FaListRecords(PSToolInterface):

    def __init__(self, args):

        super(FaListRecords, self).__init__(args)

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

        if not self.args.output in [sys.stdout, sys.stderr]:
            self.args.output = open(self.args.output, "w")


        for infile in self.args.fasta:

            for record in SeqIO.parse(infile, ftype):

                print(infile.name, record.id, len(record.seq), sep="\t", file= self.args.output)

        return True

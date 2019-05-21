import argparse

import sys

from porestat.hdf5tool import FASTQ
from porestat.tools.PTToolInterface import PSToolInterfaceFactory, PSToolInterface

from collections import OrderedDict

class FilterLengthFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(FilterLengthFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-r', '--reads', nargs='+', type=argparse.FileType('r'), help='minion read folder', required=True)

        parser.add_argument('-ub', '--max-length', type=int, help='minion read folder', required=False, default=-1)
        parser.add_argument('-lb', '--min-length', type=int, help='minion read folder', required=False, default=-1)


        def fileOpener( filename ):
            print("Opening", filename)
            open(filename, 'w').close()
            return filename

        parser.add_argument('-o', '--output', type=fileOpener, help='output location, default: std out', default=sys.stdout)
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return FilterLength(simArgs)

import os

class Environment(object):
    pass


class FilterLength(PSToolInterface):

    def __init__(self, args):

        super(FilterLength, self).__init__(args)

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

        # check https://github.com/ahcm/pyfastq_reader/blob/master/pyfastq_reader/__init__.py
        def fastq_reader_fh(infile):
            name = infile.readline().rstrip()
            while True:
                seq = ""
                for s in infile:
                    if s[0] == '+':
                        commentp = s.rstrip()
                        break
                    else:
                        seq += s.rstrip()
                qual = ""
                for q in infile:
                    if len(qual) > 0 and q[0] == '@':
                        yield name, seq, qual
                        name = q.rstrip()
                        break
                    else:
                        qual += q.rstrip()
                else:
                    yield name, seq, qual
                    return


        with open(self.args.output, 'w') as fout:
            for infile in self.args.reads:

                for rname, rseq, rqual in fastq_reader_fh(infile):

                    fq = FASTQ(rname, rseq, rqual)

                    if self.acceptFASTQ(fq):
                        fout.write(str(fq))
                        fout.write("\n")



        return True

from .PTToolInterface import PTToolInterface
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
import sys, argparse, os

class ExtractSequences(PTToolInterface):

    def __init__(self, parser, subparsers):

        super(ExtractSequences, self).__init__(parser, self.__addParser(subparsers))


    def __addParser(self, subparsers):

        parser_expls = subparsers.add_parser('seq', help='seq help')
        parser_expls.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser_expls.add_argument('-r', '--reads', type=str, help='minion read folder', required=False)
        parser_expls.add_argument('--fasta', dest='fasta', action='store_true', default=False)
        parser_expls.add_argument('--fastq', dest='fastq', action='store_true', default=True)
        parser_expls.add_argument('-e', '--experiments', nargs='+', type=str, help='run ids of experiments to be extracts', required=False)
        parser_expls.add_argument('-o', '--out', action='store', type=argparse.FileType('w'), default=sys.stdout)
        parser_expls.add_argument('-t', '--type', required=False, default =None)

        parser_expls.set_defaults(func=self.exec)

        return parser_expls

    def exec(self, args):

        folders = self.manage_folders_reads(args)
        outputstream = args.out
        endln = os.linesep

        readtype = None
        if args.type != None:
            if args.type in Fast5TYPE.str2type:
                readtype = Fast5TYPE.str2type[args.type]

        for folder in folders:

            f5folder = Fast5Directory(folder)

            for file in f5folder.collect():

                runid = file.runID()

                if args.experiments == None or len(args.experiments) == 0 or (len(args.experiments) > 0 and runid in args.experiments):

                    output = ""

                    if not args.fasta:
                        output = str(file.getFastQ(readtype))
                    else:
                        output = str(file.getFastA(readtype))

                    outputstream.write(output + endln)


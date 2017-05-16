from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory, PSToolException

from ..hdf5tool.Fast5File import Fast5Directory, Fast5TYPE
import argparse
import sys

class ExtractSequencesFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(ExtractSequencesFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser_expls = subparsers.add_parser('seq', help='seq help')
        parser_expls.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser_expls.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser_expls.add_argument('--fasta', dest='fasta', action='store_true', default=False)
        parser_expls.add_argument('--fastq', dest='fastq', action='store_true', default=True)
        parser_expls.add_argument('-e', '--experiments', nargs='+', type=str, help='run ids of experiments to be extracts', required=False)
        parser_expls.add_argument('-o', '--out', action='store', type=argparse.FileType('w'), default=sys.stdout)
        parser_expls.add_argument('-t', '--type', required=False, default =None)

        parser_expls.set_defaults(func=self._prepObj)

        return parser_expls

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return ExtractSequences(simArgs)


class Environment(object):
    pass

class ExtractSequences(ParallelPSTInterface):

    def __init__(self, args):
        super(ExtractSequences, self).__init__(args)

    def prepareEnvironment(self, args):

        oEnv = Environment()
        oEnv.out = args.out
        oEnv.experiments = args.experiments
        oEnv.fasta = args.fasta
        oEnv.fastq = args.fastq

        oEnv.readtype = None

        if args.type != None:
            if args.type in Fast5TYPE.str2type:
                oEnv.readtype = Fast5TYPE.str2type[args.type]

        return oEnv

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def execParallel(self, data, environment):

        f5folder = Fast5Directory(data)

        iFilesInFolder = 0
        collectedOutput = []

        for file in f5folder.collect():

            runid = file.runID()

            iFilesInFolder += 1

            runid = file.runID()

            if environment.experiments == None or len(environment.experiments) == 0 or (
                    len(environment.experiments) > 0 and runid in environment.experiments):

                output = None

                if not environment.fasta:
                    output = file.getFastQ(environment.readtype)
                else:
                    output = file.getFastA(environment.readtype)

                if output != None:
                    collectedOutput.append( str(output) )

        eprint("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return collectedOutput

    def joinParallel(self, existResult, newResult, oEnvironment):

        endln = os.linesep

        # TODO better endln.join(newResult) ?
        for elem in newResult:
            oEnvironment.out.write( elem + endln )

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        pass


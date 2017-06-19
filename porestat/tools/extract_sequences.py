from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory, PSToolException

from ..hdf5tool.Fast5File import Fast5Directory, Fast5TYPE
import argparse
import sys, os
from ..utils.Files import eprint


class ExtractSequencesFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(ExtractSequencesFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('seq', help='seq help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('--fasta', dest='fasta', action='store_true', default=False)
        parser.add_argument('--fastq', dest='fastq', action='store_true', default=True)
        parser.add_argument('-u', '--user_run', dest='groupByUser', action='store_true', default=False)

        parser.add_argument('-e', '--experiments', nargs='+', type=str, help='run ids of experiments to be extracted. if --user_run, give user_run_name s', required=False)
        parser.add_argument('-o', '--out', action='store', type=argparse.FileType('w'), default=sys.stdout)

        parser.add_argument('-q', '--read-type', nargs='+', dest='read_type', action=Fast5TYPEAction, default=None)

        parser.set_defaults(func=self._prepObj)

        return parser

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
        oEnv.groupByUser = args.groupByUser

        oEnv.readTypes = None
        if args.read_type != None:
            oEnv.readTypes = [ Fast5TYPE.str2type[x] for x in args.read_type]

        return oEnv

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def validFileContent(self, file, environment):
        """
        
        :param file: the Fast5File to check
        :param environment: the parallel environment / copy of args
        :return: true if the file is valid for extraction
        """


        """
        if the file does not contain reads of the desired type => continue
        """
        if not environment.readTypes is None:
            if not file.type in environment.readTypes:
                return False

        if environment.experiments == None:
            return True

        if len(environment.experiments) == 0:
            return True

        if environment.groupByUser:

            userRun = file.user_filename_input()
            return userRun in environment.experiments

        else:

            runid = file.runID()
            return runid in environment.experiments

    def execParallel(self, data, environment):

        f5folder = Fast5Directory(data)

        iFilesInFolder = 0
        collectedOutput = []

        for file in f5folder.collect():

            runid = file.runID()

            iFilesInFolder += 1

            runid = file.runID()

            if self.validFileContent(file, environment):

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

        if oEnvironment.out != sys.stdout:
            oEnvironment.out.close()


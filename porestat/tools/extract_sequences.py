from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory, PSToolException

from ..hdf5tool.Fast5File import Fast5Directory, Fast5TYPE, Fast5TYPEAction
import argparse
import sys, os
from ..utils.Files import eprint


class ExtractSequencesFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(ExtractSequencesFactory, self).__init__(parser, self._addParser(subparsers, which), which)


    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('--fasta', dest='fasta', action='store_true', default=False)
        parser.add_argument('--fastq', dest='fastq', action='store_true', default=True)
        parser.add_argument('-u', '--user_run', dest='user_run', action='store_true', default=False)

        parser.add_argument('-e', '--experiments', nargs='+', type=str, help='run ids of experiments to be extracted. if --user_run, give user_run_name s', required=False)

        def fileOpener( filename ):
            open(filename, 'w').close()
            return filename
        parser.add_argument('-o', '--output', type=fileOpener, help='output location, default: std out', default=sys.stdout)

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
        oEnv.output = args.output
        oEnv.experiments = args.experiments
        oEnv.fasta = args.fasta
        oEnv.fastq = args.fastq
        oEnv.user_run = args.user_run

        oEnv.readTypes = None
        oEnv.read_type = None

        if args.read_type != None:

            if len(args.read_type) == 1:
                oEnv.read_type = args.read_type[0]
            else:
                oEnv.read_type = None
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

        if not environment.read_type is None:
            if file.type != environment.read_type:
                return False

        if environment.experiments == None:
            return True

        if len(environment.experiments) == 0:
            return True

        if environment.user_run:

            userRun = file.user_filename_input()
            return userRun in environment.experiments

        else:

            runid = file.runID()
            return runid in environment.experiments

    def execParallel(self, datas, environment):

        for data in datas:
            f5folder = Fast5Directory(data)

            iFilesInFolder = 0
            collectedOutput = []

            for file in f5folder.collect():

                runid = file.runID()

                iFilesInFolder += 1

                if self.validFileContent(file, environment):

                    output = None

                    if not environment.fasta:
                        output = file.getFastQ(environment.read_type)
                    else:
                        output = file.getFastA(environment.read_type)

                    if output != None:
                        collectedOutput.append(str(output))

            eprint("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return collectedOutput

    def joinParallel(self, existResult, newResult, environment):

        if len(newResult) == 0:
            return existResult

        endln = os.linesep

        self.writeLinesToOutput(environment.output, [x + endln for x in newResult])
        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):
        pass

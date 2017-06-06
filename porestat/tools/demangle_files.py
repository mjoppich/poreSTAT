import os

from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory,PSToolException
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE

from ..utils.Files import makePath, fileExists, pathWritable, pathEmpty

from collections import OrderedDict, Counter
import argparse, sys, shutil


class DemangleFilesFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(DemangleFilesFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('demangle', help='expls help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)

        parser.add_argument('-q', '--read_type', nargs='+', type=str, choices=[x for x in Fast5TYPE.str2type], help='read types ('+ ",".join([x for x in Fast5TYPE.str2type]) +')')
        parser.add_argument('-u', '--user_run', dest='groupByUser', action='store_true', default=False)
        parser.add_argument('-e', '--experiments', nargs='+', type=str,
                            help='run ids of experiments to be extracted. if --user_run, give user_run_name s',
                            required=False)

        parser.add_argument('--force', action='store_true', default=False, help='Allows to choose a non-empty output folder')
        parser.add_argument('--simulate', action='store_true', default=False,
                            help='Instead of moving files, move command is printed to console')

        parser.add_argument('-c', '--chunk-size', type=int, help="how many files per subfolder. -1 for no subfolders", default=1000)

        parser.add_argument('-o', '--output', type=str, help='output location, default: std out')
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return DemangleFiles(simArgs)


class Environment(object):
    pass


class DemangleFiles(ParallelPSTInterface):

    def __init__(self, args):

        super(DemangleFiles, self).__init__(args)


    def _makePropDict(self):
        return None

    def prepareInputs(self, args):
        inputFolders =  self.manage_folders_reads(args)

        args.output = os.path.abspath(args.output)
        args.output = makePath(args.output)

        outputFolderExists = fileExists(args.output)

        if not outputFolderExists:
            os.makedirs(args.output)
            outputFolderExists = fileExists(args.output)

        outputFolderWritable = pathWritable(args.output)
        outputFolderEmpty = pathEmpty(args.output)

        if not ((outputFolderEmpty or args.force) and outputFolderExists and outputFolderWritable):
            raise PSToolException("You need to specify a writable and empty output-folder: " + str(args.output))

        return inputFolders

    def execParallel(self, data, environment):

        foundReads = []

        f5folder = Fast5Directory(data)
        iFilesInFolder = 0

        moveStatistic = Counter()
        returnMoves = {}

        for file in f5folder.collect():

            iFilesInFolder += 1

            fileExperiment = file.runID()
            fileUserRunName = file.user_filename_input()
            fileReadType = file.type

            srcPath = file.filename

            if environment.read_type != None:

                if not fileReadType in environment.read_type:
                    continue

            if environment.experiments != None:

                if not fileExperiment in environment.experiments:
                    continue

            ######### Start assembling output path
            destPath = environment.output

            if environment.groupByUser == True:
                destPath = destPath + fileUserRunName

            else:
                destPath = destPath + str(fileExperiment)

            moveStatistic[destPath] += 1

            returnMoves[ srcPath ] = destPath

        print("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return returnMoves


    def prepareEnvironment(self, args):

        env = Environment()
        env.output = args.output

        env.read_type = None
        if args.read_type != None:
            env.read_type = [Fast5TYPE[x] for x in args.read_type]

        env.groupByUser = args.groupByUser
        env.experiments = args.experiments
        env.simulate = args.simulate

        env.chunk_size = args.chunk_size

        return env


    def joinParallel(self, existResult, newResult, environment):

        if newResult == None:
            return

        if existResult is None:
            existResult = Counter()

        iChunkSize = -1
        if environment.chunk_size:
            iChunkSize = environment.chunk_size

        for source in newResult:
            destPath = newResult[source]

            if environment.chunk_size != -1:
                chunkID = existResult[destPath] % iChunkSize
                destPath = destPath + str(chunkID) + '/'

            
            self.moveFile(source, destPath, environment.simulate)

        return existResult

    def moveFile(self, src, dst, sim):

        if sim:
            print("mv " + str(src) + " " + str(dst))
        else:
            shutil.move(src, dst)


    def makeResults(self, parallelResult, oEnvironment, args):

        pass



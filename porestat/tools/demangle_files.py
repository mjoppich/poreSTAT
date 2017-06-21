import os
from enum import Enum

from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory,PSToolException
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE, Fast5TYPEAction

from ..utils.Files import makePath, fileExists, pathWritable, pathEmpty, eprint

from collections import OrderedDict, Counter
import argparse, sys, shutil

class FileHandleTypeAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, const=None, default=None, type=None, choices=None, required=False, help=None, metavar=None):
        super(FileHandleTypeAction, self).__init__(option_strings, dest, nargs, const, default, type, choices, required, help, metavar)

        self.help = 'Sets whether file gets moved or copied and must be one of {m}'.format(m=', '.join([str(x.value) for x in FileHandleType]))

    def __call__(self, parser, args, values, option_string=None):

        try:

            handleType = FileHandleType[values.upper()]
            args.__dict__[self.dest] = handleType
        except:
            raise argparse.ArgumentError(None, '{o} can not be {n}, it must be one of {m}'.format(o=str(self.option_strings), n=values, m=', '.join([str(x.name) for x in FileHandleType])))

class FileHandleType(Enum):
    MOVE='mv'
    COPY='cp'

class DemangleFilesFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(DemangleFilesFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('demangle', help='expls help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)

        parser.add_argument('-q', '--read_type', action=Fast5TYPEAction, default=None)
        parser.add_argument('-u', '--user-run', dest='groupByUser', action='store_true', default=False)
        parser.add_argument('-e', '--experiments', nargs='+', type=str,
                            help='run ids of experiments to be extracted. if --user_run, give user_run_name s',
                            required=False)

        parser.add_argument('--force', action='store_true', default=False, help='Allows to choose a non-empty output folder')
        parser.add_argument('--simulate', action='store_true', default=False,
                            help='Instead of moving files, move command is printed to console')
        parser.add_argument('--method', action=FileHandleTypeAction, default=FileHandleType.MOVE, help='')

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

        outputPath = os.path.abspath(args.output)
        args.output = makePath(outputPath)

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
                destPath = makePath(destPath + fileUserRunName)

            else:
                destPath = makePath(destPath + str(fileExperiment))

            moveStatistic[destPath] += 1

            if destPath == None:
                print(destPath)

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
        env.method = args.method

        env.chunk_size = args.chunk_size

        return env


    def joinParallel(self, existResult, newResult, environment):

        if newResult == None:
            return

        if existResult is None:
            existResult = (Counter(), Counter(), set())

        iChunkSize = -1
        if environment.chunk_size:
            iChunkSize = environment.chunk_size

        for source in newResult:
            destPath = newResult[source]
            originalDestPath = destPath

            if environment.chunk_size != -1:
                chunkID = divmod(existResult[0][destPath], iChunkSize)
                destPath = destPath + str(chunkID[0]) + '/'

            existResult[0][originalDestPath] += 1

            if not destPath in existResult[2]:
                if not environment.simulate:

                    try:
                        os.makedirs(destPath, exist_ok=True)
                        existResult[2].add(destPath)
                    except:
                        eprint(destPath + " already exists")

            existResult[1][destPath] += 1

            self.moveFile(source, destPath, environment.simulate, environment.method)

        return existResult

    def moveFile(self, src, dst, sim, method):

        if sim:
            print( method.value + " " + str(src) + " " + str(dst))
        else:

            if method == FileHandleType.MOVE:
                shutil.move(src, dst)
            elif method == FileHandleType.COPY:
                shutil.copy(src, dst)
            else:
                raise PSToolException('Invalid move file method: ' + str(method))


    def makeResults(self, parallelResult, oEnvironment, args):

        chunkFolders = sorted([x for x in parallelResult[1]])

        for x in chunkFolders:
            print(str(x) + "\t" + str(parallelResult[1][x]))


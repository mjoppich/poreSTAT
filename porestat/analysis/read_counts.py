import argparse

from ..tools.ParallelPTTInterface import ParallelPSTInterface
from ..tools.PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Parallel import Parallel as ll
from ..utils.Utils import mergeDicts
from ..utils.Stats import calcN50

import HTSeq

class ReadCountAnalysisFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(ReadCountAnalysisFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser = subparsers.add_parser('counts', help='expls help')
        parser.add_argument('-s', '--sam', nargs='+', type=str, help='alignment files')
        parser.add_argument('-g', '--gff', type=argparse.FileType("r"), help='gene annotation')
        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return ReadCountAnalysis(simArgs)



class ReadCountAnalysis(ParallelPSTInterface):

    def __init__(self, args):

        super(ReadCountAnalysis, self).__init__(args)

    def _makePropDict(self):

        props = {}

        props['COUNTS'] = {}
        props['COUNTS']['ALL'] = 0

        return props

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def execParallel(self, data, environment):

        counterRunID = {}

        f5folder = Fast5Directory(data)

        iFilesInFolder = 0

        for file in f5folder.collect():

            runid = file.runID()

            iFilesInFolder += 1

            if not runid in counterRunID:
                counterRunID[runid] = self._makePropDict()

            propDict = counterRunID[runid]

            propDict['NAMES']['USER_RUN_NAME'].add(file.user_filename_input());

            fastq = file.getFastQ()

            # if file.type == Fast5TYPE.UNKNOWN:
            #    osignal = file._get_signal()

            if fastq == None:
                propDict['TYPE'][file.type].append(0)
            else:
                propDict['TYPE'][file.type].append(len(fastq))

        print("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return counterRunID


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'FILES', 'AVG_LENGTH', 'N50', 'BASECALL_2D', 'BASECALL_1D_COMPL',
                            'BASECALL_1D', 'BASECALL_RNN_1D', 'PRE_BASECALL', 'UNKNOWN']

        if parallelResult == None:
            raise PSToolException("No Files collected")

        allobservations = {}
        for runid in parallelResult:

            allLengths = []
            countByType = Counter()

            for type in parallelResult[runid]['TYPE']:

                lengths = parallelResult[runid]['TYPE'][type]

                allLengths += lengths

                countByType[type] += len(lengths)

            fileCount = len(allLengths)
            avgLength = sum(allLengths) / fileCount
            (n50, l50) = calcN50(allLengths)
            run_user_name = parallelResult[runid]['NAMES']['USER_RUN_NAME']

            observations = {

                'RUNID': runid,
                'USER_RUN_NAME': ",".join(run_user_name),
                'FILES': fileCount,
                'AVG_LENGTH': avgLength,
                'N50': n50,
                'BASECALL_2D': countByType[self.str2fileType['BASECALL_2D']],
                'BASECALL_1D_COMPL': countByType[self.str2fileType['BASECALL_1D_COMPL']],
                'BASECALL_1D': countByType[self.str2fileType['BASECALL_1D']],
                'BASECALL_RNN_1D': countByType[self.str2fileType['BASECALL_RNN_1D']],
                'PRE_BASECALL': countByType[self.str2fileType['PRE_BASECALL']],
                'UNKNOWN': countByType[self.str2fileType['UNKNOWN']]
            }

            allobservations[runid] = observations


        sortedruns = sorted([x for x in allobservations])

        print("\t".join(makeObservations))

        for runid in sortedruns:

            allobs = []
            for x in makeObservations:
                allobs.append(str(allobservations[runid][x]))

            print("\t".join(allobs))



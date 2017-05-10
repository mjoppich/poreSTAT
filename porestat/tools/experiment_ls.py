from .ParallelPTTInterface import ParallelPTTInterface
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Parallel import Parallel as ll
from ..utils.Utils import mergeDicts


class Experiment_ls(ParallelPTTInterface):

    def __init__(self, parser, subparsers):

        super(Experiment_ls, self).__init__(parser, self.__addParser(subparsers))

        self.fileType2Str = {

            Fast5TYPE.BASECALL_2D: 'BASECALL_2D',
            Fast5TYPE.BASECALL_1D_COMPL: 'BASECALL_1D_COMPL',
            Fast5TYPE.BASECALL_1D: 'BASECALL_1D',
            Fast5TYPE.BASECALL_RNN_1D: 'BASECALL_RNN_1D',
            Fast5TYPE.PRE_BASECALL: 'PRE_BASECALL',
            Fast5TYPE.UNKNOWN: 'UNKNOWN'
         }

        self.str2fileType = {self.fileType2Str[x]: x for x in self.fileType2Str}

    def __addParser(self, subparsers):

        parser_expls = subparsers.add_parser('expls', help='expls help')
        parser_expls.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser_expls.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser_expls.set_defaults(func=self.exec)

        return parser_expls

    def _makePropDict(self):

        props = {}

        props['TYPE'] = {}

        for ftype in self.fileType2Str:
            props['TYPE'][ftype] = []

        props['NAMES'] = {}
        props['NAMES']['USER_RUN_NAME'] = set()

        return props

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def execParallel(self, procID, environment, data):

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
            n50 = self._calcN50(allLengths)
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

    def _calcN50(self, lengths):

        slengths = sorted(lengths, reverse=True)
        totalLength = sum(slengths)

        currentLength = 0
        neededLength = totalLength / 2.0

        print(neededLength)

        n50Value = 0
        L50 = 0

        for i in range(0, len(slengths)):

            x = slengths[i]

            currentLength += x
            L50 += 1

            if currentLength >= neededLength:
                n50Value = x

                print("N50 value is " + str(x) + " for length " + str(totalLength) + " (in half: " + str(neededLength) + " ) and L50: " + str(i) + " " + str(len(slengths)))
                break


        return n50Value

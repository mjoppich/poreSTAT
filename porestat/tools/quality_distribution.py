from .ParallelPTTInterface import ParallelPTTInterface
from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Parallel import Parallel as ll
from ..utils.Utils import mergeDicts, mergeCounter


class QualityDistribution(ParallelPTTInterface):

    def __init__(self, parser, subparsers):

        super(QualityDistribution, self).__init__(parser, self.__addParser(subparsers))

        self.qualTypes = [chr(x) for x in range(ord('!'), ord('~')+1)]

    def __addParser(self, subparsers):

        parser_expls = subparsers.add_parser('qual_dist', help='expls help')
        parser_expls.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser_expls.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser_expls.add_argument('-p', '--plot', nargs='?', type=bool, const=True, default=False, help='issue plot?', required=False)
        parser_expls.set_defaults(func=self.exec)

        return parser_expls

    def _makePropDict(self):

        propDict = {}
        propDict['QUALS'] = Counter()
        propDict['USER_RUN_NAME'] = set()
        propDict['READ_COUNT'] = 0

        return propDict

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
            propDict['READ_COUNT'] += 1
            propDict['USER_RUN_NAME'].add( file.user_filename_input() )

            fastq = file.getFastQ()

            if fastq != None:

                for x in fastq.qual:
                    propDict['QUALS'][x] += 1


        print("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return counterRunID


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'FILES', 'TOTAL_BASES']
        for x in self.qualTypes:
            makeObservations.append(x)

        for x in self.qualTypes:
            makeObservations.append(x + "%")

        allobservations = {}
        for runid in parallelResult:

            props = parallelResult[runid]

            run_user_name = props['USER_RUN_NAME']
            fileCount = props['READ_COUNT']

            nuclCounts = props['QUALS']

            observations = {

                'RUNID': runid,
                'USER_RUN_NAME': ",".join(run_user_name),
                'FILES': fileCount
            }

            allNucl = 0
            for x in self.qualTypes:
                observations[x] = nuclCounts[x]
                allNucl += nuclCounts[x]

            observations['TOTAL_BASES'] = allNucl

            for x in self.qualTypes:
                observations[x + "%"] = nuclCounts[x] / allNucl

            allobservations[runid] = observations

        sortedruns = sorted([x for x in allobservations])

        print("\t".join(makeObservations))

        for runid in sortedruns:

            allobs = []
            for x in makeObservations:
                allobs.append(str(allobservations[runid][x]))

            print("\t".join(allobs))


from .ParallelPTTInterface import ParallelPSTInterface
from .PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

import argparse
class QualityPositionFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(QualityPositionFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):

        parser_expls = subparsers.add_parser('qual_pos', help='expls help')
        parser_expls.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser_expls.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser_expls.add_argument('-p', '--plot', '--out', action='store', type=argparse.FileType('w'), default=None)
        parser_expls.set_defaults(func=self._prepObj)

        return parser_expls

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return QualityPosition(simArgs)

class QualityPosition(ParallelPSTInterface):

    def __init__(self, args):

        super(QualityPosition, self).__init__( args )

        self.qualTypes = [chr(x) for x in range(ord('!'), ord('~')+1)]



    def _makePropDict(self):

        propDict = {}
        propDict['QUALS'] = {}
        propDict['USER_RUN_NAME'] = set()
        propDict['READ_COUNT'] = 0

        return propDict

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def execParallel(self, procID, environment, data):

        counterRunID = {}

        f5folder = Fast5Directory(data)

        iFilesInFolder = 0

        print("Folder started: " + f5folder.path)

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

                qualDict = propDict['QUALS']

                for i in range(0, len(fastq.qual)):

                    if not i in qualDict:
                        qualDict[i] = Counter()

                    qualDict[i][fastq.qual[i]] += 1

                propDict['QUALS'] = qualDict


        print("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return counterRunID


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'FILES']
        for x in self.qualTypes:
            makeObservations.append(x)

        for x in self.qualTypes:
            makeObservations.append(x + "%")

        allobservations = {}
        for runid in parallelResult:

            props = parallelResult[runid]

            run_user_name = props['USER_RUN_NAME']
            fileCount = props['READ_COUNT']

            qualCounter = props['QUALS']

            observations = {

                'RUNID': runid,
                'USER_RUN_NAME': ",".join(run_user_name),
                'FILES': fileCount,
                'QUALPOS': qualCounter
            }

            allobservations[runid] = observations

        sortedruns = sorted([x for x in allobservations])

        print("\t".join(makeObservations))

        for runid in sortedruns:

            allobs = []
            for x in makeObservations:
                allobs.append(str(allobservations[runid][x]))

            print("\t".join(allobs))

            # make plot for runid

            # pos -> qual -> count
            qualCounter = observations['QUALPOS']

            self.plotQualCounter(qualCounter)


    def plotQualCounter(self, qualCounter):


        foundLengths = set()
        for length in qualCounter:
            foundLengths.add(length)

        foundLengths = sorted([foundLengths])
        maxLength = max(foundLengths)
        step = maxLength / 10;

        fig, axes = plt.subplots(nrows=1, ncols=10, figsize=(6, 6))

        for i in range(0, 10):

            stepMin = i*step
            stepMax = stepMin + step

            allData = Counter()
            for pos in qualCounter:

                if stepMin <= pos and pos < stepMax:

                    for k in qualCounter[pos]:
                        allData[k] += qualCounter[pos][k]

            dataVal = []
            dataPos = []

            for x in allData:
                dataVal.append(allData[x])
                dataPos.append(x)

            axes[0, i].violinplot(dataVal, dataPos, points=20, widths=0.3, showmeans=True, showextrema=True, showmedians=True)





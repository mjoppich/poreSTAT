from porestat.plots.plotconfig import PlotConfig

from .ParallelPTTInterface import ParallelPSTReportableInterface
from .PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

class NucleotideDistributionFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(NucleotideDistributionFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('nuc_dist', help='expls help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        simArgs = self._makeArguments(args)

        return NucleotideDistribution(simArgs)

class NucleotideDistribution(ParallelPSTReportableInterface):

    def __init__(self, args):

        super(NucleotideDistribution, self).__init__( args )

        self.nucTypes = [
            'A','C','T','G','N'
         ]



    def _makePropDict(self):

        propDict = {}
        propDict['NUCS'] = Counter()
        propDict['USER_RUN_NAME'] = set()
        propDict['READ_COUNT'] = 0

        return propDict

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def handleEntity(self, fileObj, localEnv, globalEnv):

        runid = fileObj.runID()

        if not runid in localEnv:
            localEnv[runid] = self._makePropDict()

        propDict = localEnv[runid]
        propDict['READ_COUNT'] += 1
        propDict['USER_RUN_NAME'].add(fileObj.user_filename_input())

        fastq = fileObj.getFastQ()

        if fastq != None:

            for x in fastq.seq:
                propDict['NUCS'][x] += 1

        return localEnv


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'FILES', 'TOTAL_BASES']

        for x in self.nucTypes:
            makeObservations.append(x)

        for x in self.nucTypes:
            makeObservations.append(x + "%")

        allobservations = {}
        for runid in parallelResult:

            props = parallelResult[runid]

            run_user_name = props['USER_RUN_NAME']
            fileCount = props['READ_COUNT']

            nuclCounts = props['NUCS']

            observations = {

                'RUNID': runid,
                'USER_RUN_NAME': ",".join(run_user_name),
                'FILES': fileCount
            }

            allNucl = 0
            for x in self.nucTypes:
                observations[x] = nuclCounts[x]
                allNucl += nuclCounts[x]

            observations['TOTAL_BASES'] = allNucl

            for x in self.nucTypes:
                observations[x + "%"] = nuclCounts[x] / allNucl

            allobservations[runid] = observations

        sortedruns = sorted([x for x in allobservations])

        print("\t".join(makeObservations))

        for runid in sortedruns:

            allobs = []
            for x in makeObservations:
                allobs.append(str(allobservations[runid][x]))

            print("\t".join(allobs))

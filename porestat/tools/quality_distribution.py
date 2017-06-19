from collections import OrderedDict

from ..plots.plotconfig import PlotConfig
from ..plots.poreplot import PorePlot

from .ParallelPTTInterface import ParallelPSTReportableInterface
from .PTToolInterface import PSToolInterfaceFactory

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

class QualityDistributionFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(QualityDistributionFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('qual_dist', help='expls help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-np', '--no-plot', nargs='?', type=bool, const=True, default=False, help='set if no plot should be issued', required=False)
        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        simArgs = self._makeArguments(args)

        return QualityDistribution(simArgs)

class QualityDistribution(ParallelPSTReportableInterface):

    def __init__(self, args):

        super(QualityDistribution, self).__init__( args )

        self.qualTypes = [chr(x) for x in range(ord('!'), ord('~')+1)]


    def _makePropDict(self):

        propDict = {}
        propDict['QUALS'] = Counter()
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

            for x in fastq.qual:
                propDict['QUALS'][x] += 1

        return localEnv


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        statObservations = ['RUNID', 'USER_RUN_NAME', 'FILES', 'TOTAL_BASES']
        absQualities = [x for x in self.qualTypes]
        relQualitites = [str(x) + "%" for x in self.qualTypes]

        makeObservations = statObservations + absQualities + relQualitites

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
            for x in absQualities:
                observations[x] = nuclCounts[x]
                allNucl += nuclCounts[x]

            observations['TOTAL_BASES'] = allNucl

            for x in absQualities:
                observations[x+'%'] = nuclCounts[x] / allNucl

            key = ",".join(run_user_name) if self.hasArgument('groupByRunName', args) and args.groupByRunName else runid

            if key in allobservations:
                allobservations[key] = mergeDicts(allobservations[key], observations)
            else:
                allobservations[key] = observations


        sortedruns = sorted([x for x in allobservations])

        print("\t".join(makeObservations))

        absPlotData = {}
        relPlotData = {}

        for runid in sortedruns:

            allobs = []
            absCounts = OrderedDict()
            relCounts = OrderedDict()

            for x in statObservations:
                allobs.append(str(allobservations[runid][x]))

            for x in absQualities:
                allobs.append(str(allobservations[runid][x]))

                absCounts[x] = allobservations[runid][x]
                relCounts[x] = allobservations[runid][str(x)+'%']

            absPlotData[runid] = absCounts
            relPlotData[runid] = relCounts

            print("\t".join(allobs))
        
        PorePlot.plotBars(absPlotData, "Quality Distribution (abs. values)", "Qualities", "Count (nucl.)", pltcfg=args.pltcfg)
        PorePlot.plotBars(relPlotData, "Quality Distribution (perc)", "Qualities", "percentage", pltcfg=args.pltcfg)

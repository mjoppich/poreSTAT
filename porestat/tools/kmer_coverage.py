import argparse
import math
import sys
from collections import defaultdict, OrderedDict

from Bio import SeqIO
from porestat.utils.DataFrame import DataFrame, DataRow

from ..plots.plotconfig import PlotConfig
from ..plots.poreplot import PorePlot, PlotDirectionTYPE

from .ParallelPTTInterface import ParallelPSTReportableInterface
from .PTToolInterface import PSToolInterfaceFactory, PSToolException
from ..utils.Stats import calcN50

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Utils import mergeDicts, mergeCounter

class KmerHistogramFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers, which):

        super(KmerHistogramFactory, self).__init__(parser, self._addParser(subparsers, which), which)



    def _addParser(self, subparsers, which):
        parser = subparsers.add_parser(which, help=which+' help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)
        parser.add_argument('-p', '--plot', nargs='?', type=bool, const=True, default=False, help='issue plot?', required=False)
        parser.add_argument('-u', '--user-run', dest='user_run', action='store_true', default=False)
        parser.add_argument('-q', '--read-type', dest='read_type', action='store_true', default=False, help='add type subplots')
        parser.add_argument('-k', '--kmer', dest='k', type=int, default=5, help='add type subplots')
        parser.add_argument('-mc', '--mostcommon', dest='mc', type=int, default=10, help='add type subplots')

        parser.add_argument('--reference', dest='reference', type=argparse.FileType('r'), help='fasta reference for kmer distribution', default=None)


        parser.add_argument('-c', '--combined', dest='combineRuns', action='store_true', default=False)
        parser.add_argument('-v', '--violin', dest='violin', action='store_true', default=False)

        parser = PlotConfig.addParserArgs(parser)

        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):
        simArgs = self._makeArguments(args)
        simArgs.pltcfg = PlotConfig.fromParserArgs(simArgs)

        return KmerHistogram(simArgs)

class KmerHistogram(ParallelPSTReportableInterface):

    def __init__(self, args):

        super(KmerHistogram, self).__init__( args )

        if args.k > 10:
            sys.stderr.write('WARNING: You specified k={k} which can result in a huge memory consumption. Consider lowering your k if this run crashes.'.format(k=args.k))
            sys.stderr.write(
                'Try tsxCount ( http://github.com/mjoppich/tsxcount ) for fast and efficient k-mer counting.')


    def _makePropDict(self):

        propDict = {}
        propDict['USER_RUN_NAME'] = set()
        propDict['READ_COUNT'] = 0
        propDict['KMERS'] = Counter()

        return propDict

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)


    @classmethod
    def calcKmers(cls, seqStr, k):

        retCounter = Counter()

        if seqStr != None:

            for i in range(0, len(seqStr)-k):
                kmer = seqStr[i:i+k]
                assert(len(kmer)==k)

                retCounter[kmer] += 1

        return retCounter


    def handleEntity(self, fileObj, localEnv, globalEnv):

        runid = fileObj.runID()

        if not runid in localEnv:
            localEnv[runid] = self._makePropDict()

        propDict = localEnv[runid]
        propDict['READ_COUNT'] += 1
        propDict['USER_RUN_NAME'].add(fileObj.user_filename_input())

        fastq = fileObj.getFastQ()

        if fastq != None:
            propDict['KMERS'] += self.calcKmers(fastq.seq, self.args.k)

        return localEnv

    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

        if parallelResult == None:
            raise PSToolException('No valid result generated.')

        makeObservations = ['RUNID', 'USER_RUN_NAME', 'KMER', 'COUNT', 'COUNT_MEAN', 'COUNT_RMSD']

        """

        ONLY IF REFERENCE IS GIVEN => calc reference distribution

        """

        if args.reference != None:

            refKmerCounts = Counter()
            for record in SeqIO.parse(args.reference, "fasta"):
                refKmerCounts += self.calcKmers(str(record.seq), self.args.k)

            parallelResult['fasta_reference'] = {
                'USER_RUN_NAME': ['fasta_reference'],
                'READ_COUNT': 1,
                'KMERS': refKmerCounts
            }

        allobservations = {}
        for runid in parallelResult:

            props = parallelResult[runid]

            run_user_name = props['USER_RUN_NAME']
            fileCount = props['READ_COUNT']
            kmerObservations = props['KMERS']

            observations = {

                'RUNID': runid,
                'USER_RUN_NAME': ",".join(run_user_name),
                'FILES': fileCount,
                'KMERCOUNTS': kmerObservations
            }

            key = self.makeKey(run_user_name, args, runid)

            if key in allobservations:
                allobservations[key] = mergeDicts(allobservations[key], observations)
            else:
                allobservations[key] = observations

        sortedruns = sorted([x for x in allobservations])

        plotData = {}

        if self.hasArgument('combineRuns', args) and args.combineRuns:

            labels = []

            kmerData = Counter()

            for runid in sortedruns:
                kmerData += allobservations[runid]['KMERCOUNTS']
                labels.append(runid)


                args.pltcfg.makeTable(self.dfSummary(allobservations[runid], self.args.mc))

            plotLabel = ",".join(labels)
            plotData[ plotLabel ] = [kmerData[x] for x in kmerData]

        else:

            for runid in sortedruns:

                args.pltcfg.makeTable(self.dfSummary(allobservations[runid], self.args.mc))
                kmerData = allobservations[runid]['KMERCOUNTS']

                plotData[runid] =  [kmerData[x] for x in kmerData]

        if self.hasArgument('read_type', args) and args.read_type:

            newPlotData = {}

            for runid in plotData:
                lengthsByType = defaultdict(list)

                for x in plotData[runid]:
                    lengthsByType[x[1]].append(x[0])

                for readtype in lengthsByType:
                    newid = runid + "_" + readtype
                    newPlotData[newid] = lengthsByType[readtype]

            plotData = newPlotData




        if self.hasArgument('violin', args) and args.violin:
            PorePlot.plotViolin(plotData, None, 'Kmer Histogram', pltcfg=args.pltcfg, plotDirection=PlotDirectionTYPE.HORIZONTAL)
        else:
            PorePlot.plotHistogram(plotData, None, 'Kmer Histogram for ', xlabel="k-mer Frequency", ylabel="k-mer count", pltcfg=args.pltcfg, plotDirection=PlotDirectionTYPE.VERTICAL)
            PorePlot.plotCumHistogram(plotData, None, 'Kmer Histogram for ', xlabel="k-mer Frequency", bins=-1,
                                   ylabel="Number k-mers", pltcfg=args.pltcfg, xLogAxis=False, normed=False, plotDirection=PlotDirectionTYPE.VERTICAL)



    @classmethod
    def summarizeKmers(cls, thisObservation, mc=10):

        # ['RUNID', 'USER_RUN_NAME', 'KMER', 'COUNT']
        obsInfo = {
            'RUNID': thisObservation['RUNID'],
            'USER_RUN_NAME': thisObservation['USER_RUN_NAME'],
        }

        kmers = thisObservation['KMERCOUNTS']

        allobs = []

        for kmer, count in kmers.most_common(mc):

            kmercounts = [kmers[x] for x in kmers]
            mean = sum(kmercounts) / len(kmers)
            var = math.sqrt(sum([(xi - mean) ** 2 for xi in kmercounts]) / (len(kmercounts) - 1))

            thisobs = OrderedDict([
                ('RUNID', obsInfo['RUNID']),
                ('RUN_NAME', obsInfo['USER_RUN_NAME']),
                ('KMER', str(kmer)),
                ('KMER_COUNT', str(count)),
                ('KMER_COUNT_MEAN', str(mean)),
                ('KMER_COUNT_RMSD', str(var))

            ])

            allobs.append(thisobs)

        return allobs

    @classmethod
    def dfSummary(cls, thisObservation, mc=10):

        df = DataFrame()

        madeObs = cls.summarizeKmers(thisObservation, mc)

        for idx, obs in enumerate(madeObs):

            if idx == 0: # header
                cols = [x for x in obs]

                df.addColumns(cols)

            dfrow = DataRow.fromDict(obs)

            df.addRow(dfrow)

        return df




    @classmethod
    def printObservation(cls, thisObservation, mc=10):

        madeObs = cls.summarizeKmers(thisObservation, mc)
        allobs = []

        for idx, obs in enumerate(madeObs):

            if idx == 0: # header
                allobs.append("\t".join([x for x in obs]))

            allobs.append("\t".join([obs[x] for x in obs]))

        print("\n".join(allobs))


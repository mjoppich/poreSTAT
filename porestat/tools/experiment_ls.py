import argparse
from collections import OrderedDict

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

from .ParallelPTTInterface import ParallelPSTReportableInterface
from .PTToolInterface import PSToolInterfaceFactory,PSToolException

from ..hdf5tool.Fast5File import Fast5File, Fast5Directory, Fast5TYPE
from collections import Counter
from ..utils.Parallel import Parallel as ll
from ..utils.Utils import mergeDicts
from ..utils.Stats import calcN50

class ExperimentLsFactory(PSToolInterfaceFactory):

    def __init__(self, parser, subparsers):

        super(ExperimentLsFactory, self).__init__(parser, self._addParser(subparsers))


    def _addParser(self, subparsers):
        parser = subparsers.add_parser('expls', help='expls help')
        parser.add_argument('-f', '--folders', nargs='+', type=str, help='folders to scan', required=False)
        parser.add_argument('-r', '--reads', nargs='+', type=str, help='minion read folder', required=False)

        parser.add_argument('-x', '--xlsx', type=str, default=None, help='Save output to excel file (xlsx)')


        parser.set_defaults(func=self._prepObj)

        return parser

    def _prepObj(self, args):

        simArgs = self._makeArguments(args)

        return ExperimentLs(simArgs)



class ExperimentLs(ParallelPSTReportableInterface):

    def __init__(self, args):

        super(ExperimentLs, self).__init__(args)

        self.fileTypes = OrderedDict([(x, x.value) for x in Fast5TYPE])

        #self.str2fileType = {self.fileType2Str[x]: x for x in self.fileType2Str}

    def _makePropDict(self):

        props = {}

        props['TYPE'] = {}

        for ftype in self.fileTypes:
            props['TYPE'][ftype] = []

        props['NAMES'] = {}
        props['NAMES']['USER_RUN_NAME'] = set()

        return props

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def handleEntity(self, fileObj, localEnv, globalEnv):

        runid = fileObj.runID()

        if not runid in localEnv:
            localEnv[runid] = self._makePropDict()

        propDict = localEnv[runid]

        propDict['NAMES']['USER_RUN_NAME'].add(fileObj.user_filename_input());

        fastq = fileObj.getFastQ()

        # if file.type == Fast5TYPE.UNKNOWN:
        #    osignal = file._get_signal()

        if fastq == None:
            propDict['TYPE'][fileObj.type].append(0)
        else:
            propDict['TYPE'][fileObj.type].append(len(fastq))

        return localEnv


    def joinParallel(self, existResult, newResult, oEnvironment):

        if existResult == None:
            existResult = {}

        existResult = mergeDicts(existResult, newResult)

        return existResult


    def makeResults(self, parallelResult, oEnvironment, args):

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

            observations = OrderedDict([
                ('RUNID', runid),
                ('USER_RUN_NAME', ",".join(run_user_name)),
                ('FILES', fileCount),
                ('AVG_LENGTH', avgLength),
                ('N50', n50),
            ])
            observations = mergeDicts(observations, OrderedDict( [ (s, countByType[ t ]) for (t, s) in self.fileTypes.items() ] ), OrderedDict)

            allobservations[runid] = observations


        if len(allobservations) > 0:

            exobs = list(allobservations.items())[0][1]
            headerNames = [x for x in exobs]

            sortedruns = sorted([x for x in allobservations])

            print("\t".join(headerNames))

            for runid in sortedruns:

                allobs = []
                for x in allobservations[runid]:
                    allobs.append(str(allobservations[runid][x]))

                print("\t".join(allobs))

            if args.xlsx:
                outFrame = DataFrame()
                outFrame.addColumns(headerNames)

                for runid in sortedruns:
                    row = DataRow.fromDict( allobservations[runid] )
                    outFrame.addRow( row )

                outFrame.export(args.xlsx, ExportTYPE.XLSX)


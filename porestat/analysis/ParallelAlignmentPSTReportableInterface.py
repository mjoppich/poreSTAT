import time

import pysam

import HTSeq

from porestat.utils.Parallel import MapReduce
from porestat.utils import eprint
from ..tools.ParallelPTTInterface import ParallelPSTInterface


class ParallelAlignmentPSTReportableInterface(ParallelPSTInterface):
    def __init__(self, args):

        super(ParallelAlignmentPSTReportableInterface, self).__init__(args)

    def _createLocalEnvironment(self):
        return {}

    def handleEntity(self, headers, samStrings, envs):
        pass


    def exec(self):

        iStart = time.time()
        inputs = self.prepareInputs(self.args)
        iEnd = time.time()
        eprint("Preparing Inputs: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd - iStart))))

        print(inputs)
        environment = self.prepareEnvironment(self.args)
        ll = MapReduce(8)
        llResults = ll.exec( inputs, self.execParallel, environment, self.chunkSize, self.joinParallel)
        iEnd = time.time()
        eprint("Exec Parallel: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd - iStart))))

        self.makeResults(llResults, None, self.args)
        iEnd = time.time()
        eprint("Results Time: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd-iStart))))

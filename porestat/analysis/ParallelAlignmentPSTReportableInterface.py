import time

import pysam

import HTSeq

from porestat.analysis.ParallelAlignment import MapReduceAlignment
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

        llResults = self.execParallel(inputs, None)
        iEnd = time.time()
        eprint("Exec Parallel: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd - iStart))))

        self.makeResults(llResults, None, self.args)
        iEnd = time.time()
        eprint("Results Time: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd-iStart))))

from .PTToolInterface import PSToolInterface
from ..utils.Parallel import MapReduce
from ..utils.Utils import eprint
import time

class ParallelPSTInterface(PSToolInterface):

    def __init__(self, args):

        super(ParallelPSTInterface, self).__init__(args)

        self.chunkSize = 1


    def exec(self):

        iStart = time.time()
        environment = self.prepareEnvironment(self.args)
        inputs = self.prepareInputs(self.args)

        ll = MapReduce(4)

        result = ll.exec( inputs, self.execParallel, environment, self.chunkSize, self.joinParallel)

        self.makeResults(result, environment, self.args)
        iEnd = time.time()

        eprint("Execution Time: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd-iStart))))

    def prepareEnvironment(self, args):
        return self._makeArguments(args)

    def prepareInputs(self, args):
        return self.manage_folders_reads(args)

    def execParallel(self, data, environment):
        pass

    def joinParallel(self, existResult, newResult, oEnvironment):
        pass

    def makeResults(self, parallelResult, environment, args):
        pass



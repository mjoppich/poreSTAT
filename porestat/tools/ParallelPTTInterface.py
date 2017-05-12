from .PTToolInterface import PTToolInterface
from ..utils.Parallel import Parallel as ll
from ..utils.Utils import eprint
import time

class ParallelPTTInterface(PTToolInterface):

    def __init__(self, parser, subparsers):

        super(ParallelPTTInterface, self).__init__(parser, subparsers)

        self.chunkSize = 1


    def exec(self, args):

        iStart = time.time()
        environment = self.prepareEnvironment(args)
        inputs = self.prepareInputs(args)

        result = ll.mapReduce(4, inputs, self.execParallel, environment, self.chunkSize, self.joinParallel)

        self.makeResults(result, environment, args)
        iEnd = time.time()

        eprint("Execution Time: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd-iStart))))

    def prepareEnvironment(self, args):
        return None

    def prepareInputs(self, args):
        pass

    def execParallel(self, procID, environment, data):
        pass

    def joinParallel(self, existResult, newResult, oEnvironment):
        pass

    def makeResults(self, parallelResult, environment, args):
        pass



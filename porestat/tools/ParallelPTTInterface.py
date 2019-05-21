from porestat.hdf5tool import Fast5Directory
from .PTToolInterface import PSToolInterface
from ..utils.Parallel import MapReduce
from porestat.utils import eprint
import time


class ParallelPSTInterface(PSToolInterface):

    def __init__(self, args):

        super(ParallelPSTInterface, self).__init__(args)

        self.chunkSize = 1


    def exec(self):

        iStart = time.time()
        inputs = self.prepareInputs(self.args)
        environment = self.prepareEnvironment(self.args)

        ll = MapReduce(4)
        result = ll.exec( inputs, self.execParallel, environment, self.chunkSize, self.joinParallel)

        iEnd = time.time()
        eprint("Making Results: " + str(time.strftime('%H:%M:%S [HH:MM:SS]', time.gmtime(iEnd - iStart))))

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

    def writeLinesToOutput(self, outFile, lines, mode='a'):

        if type(outFile) == str:
            file = open(outFile, mode)
            file.writelines(lines)
            file.flush()
            file.close()

        else:
            for line in lines:
                outFile.write(line)

            outFile.flush()

    def closeOutput(self, outFile):

        if type(outFile) != str:
            outFile.close()


class ParallelPSTReportableInterface(ParallelPSTInterface):

    def __init__(self, args):

        super(ParallelPSTReportableInterface, self).__init__(args)

    def _createLocalEnvironment(self):
        return {}

    def handleEntity(self, fileObj, localEnv, globalEnv):
        pass

    def execParallel(self, data, environment):
        iFilesInFolder = 0

        localEnv = self._createLocalEnvironment()

        for folder in data:

            f5folder = Fast5Directory(folder)

            for file in f5folder.collect():
                localEnv = self.handleEntity(file, localEnv, environment)
                iFilesInFolder += 1

            print("Folder done: " + f5folder.path + " [Files: " + str(iFilesInFolder) + "]")

        return localEnv



import HTSeq

from ..tools.ParallelPTTInterface import ParallelPSTInterface


class ParallelAlignmentPSTReportableInterface(ParallelPSTInterface):
    def __init__(self, args):

        super(ParallelAlignmentPSTReportableInterface, self).__init__(args)

    def _createLocalEnvironment(self):
        return {}

    def handleEntity(self, fileObj, localEnv, globalEnv):
        pass

    def execParallel(self, data, environment):

        retData = []
        for alignFile in data:

            if alignFile.name.endswith(".bam"):
                opener = HTSeq.BAM_Reader
            else:
                opener = HTSeq.SAM_Reader

            iProcessedAlignments = 0
            localEnv = self._createLocalEnvironment()

            for readAlignment in opener(alignFile):
                localEnv = self.handleEntity(readAlignment, localEnv, environment)
                iProcessedAlignments += 1

            print("In File " + str(alignFile.name) + " " + str(iProcessedAlignments) + " have been processed")

            retData.append( (alignFile.name, localEnv) )

        return retData
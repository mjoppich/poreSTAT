import time
from pathos import multiprocessing as mp
import itertools

class MapReduceAlignment:

    def __init__(self, procs = 4):
        self.pool = None
        self.nprocs = procs

    def exec(self, samFile, oFunc, sEnvironment, chunkSize = 1000, pReduceFunc = None):

        self.pool = mp.ProcessPool(self.nprocs)
        allResults = []

        resultObj = None


        samHeader = samFile.header

        hasEnded = False

        while len(allResults) < self.nprocs:

            samChunk = []
            for aln in samFile:

                samChunk.append(aln.to_string())

                if len(samChunk) >= chunkSize:
                    break


            if len(samChunk) != chunkSize:
                hasEnded = True

            allResults.append( self.pool.apipe( oFunc, samHeader, samChunk, sEnvironment ) )


        if hasEnded:
            self.pool.close()


        while len(allResults) > 0:

            i=0
            while i < len(allResults):

                if allResults[i].ready():

                    result = allResults[i].get()

                    if pReduceFunc != None:

                        resultObj = pReduceFunc(resultObj, result, sEnvironment)

                    else:

                        if resultObj == None:
                            resultObj = []

                        resultObj.append(result)

                    del allResults[i]



                    """
                    FILL UP THE QUEUE
                    """

                    if not hasEnded:
                        while len(allResults) < self.nprocs:

                            samChunk = []
                            for aln in samFile:

                                samChunk.append(aln.to_string())

                                if len(samChunk) >= chunkSize:
                                    break

                            if len(samChunk) != chunkSize:
                                hasEnded = True

                        if hasEnded:
                            self.pool.close()

                else:
                    i += 1

            time.sleep(0.5)

        print("Pool Join")
        self.pool.join()
        print("Pool Clear")
        self.pool.clear()
        print("Pool Closed")

        return resultObj


    @classmethod
    def chunkIterable(cls, iterable, size):

        it = iter(iterable)
        item = list(itertools.islice(it, size))

        while item:
            yield item
            item = list(itertools.islice(it, size))
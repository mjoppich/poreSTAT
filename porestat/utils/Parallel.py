__author__ = 'joppich'

import queue
import threading
import time
import sys
import os
import multiprocessing
import pickle

from .Utils import eprint

class Parallel:

    datamanager = None

    class Thread (threading.Thread):
        def __init__(self, threadID, oQueue, oLock, oFunc):
            threading.Thread.__init__(self)
            self.threadID = threadID
            self.oQueue = oQueue
            self.oLock = oLock

            self.process = oFunc

            self.bFinished = False

        def run(self):
            eprint ( "Starting Thread" + str(self.threadID))
            self.process_data()
            eprint ("Exiting Thread" + str(self.threadID))

        def process_data(self):

            while not self.bFinished:
                self.oLock.acquire()
                if not self.oQueue.empty():
                    data = self.oQueue.get()
                    self.oLock.release()

                    self.process(self.threadID, data)

                else:
                    self.oLock.release()

                time.sleep(0.1)

        def cancel(self):
            self.bFinished = True

        def reset(self):
            self.bFinished = False

    @classmethod
    def foreach(cls, iThreads, oIterable, oFunc):

        """

        Thread based parallelism

        :param iThreads:
        :param oIterable:
        :param oFunc:
        :return:
        """

        queueLock = threading.Lock()
        workQueue = queue.Queue( len(oIterable) )
        vThreads = []

        bFinished = False

        # Create new threads
        for iThreadID in range(0, iThreads):

            thread = Parallel.Thread(iThreadID, workQueue,queueLock, oFunc)

            thread.start()
            vThreads.append(thread)

        # now fill the queue
        queueLock.acquire()
        for oElem in oIterable:
            workQueue.put( oElem )
        queueLock.release()

        # wait for all elements processed # TODO why busy waiting?
        while not workQueue.empty():
            pass

        # finish all threads
        for oThread in vThreads:
            oThread.cancel()
            oThread.join()

    class MyProcess(multiprocessing.Process):

        def __init__(self, threadID, vChunk, vReturnQueue, oReturnQueueLock, sEnvironment, oFunc):
            super(Parallel.MyProcess, self).__init__()
            self.procID = threadID
            self.vChunk = vChunk
            self.process = oFunc
            self.vReturnQueue = vReturnQueue
            self.oReturnQueueLock = oReturnQueueLock

            self.bFinished = False
            self.sEnvironment = sEnvironment

        def __str__(self):

            return "MyProcess " + str(self.procID) + " " + str(self.vChunk)

        def run(self):

            #eprint("Running " + str(self))

            for data in self.vChunk:

                oReturn = self.process(self.procID, self.sEnvironment, data)
                self.vReturnQueue.put(oReturn)
                self.vReturnQueue.join()

            #eprint("Closing " + str(self))
            self.vReturnQueue.close()
            self.vReturnQueue.join_thread()

            #eprint("Joined " + str(self))

            self.bFinished = True

            iCode = 0
            #os.kill(os.getpid(), 9)
            return iCode

    @classmethod
    def getDataManager(cls):

        if cls.datamanager == None:
            cls.datamanager = multiprocessing.Manager()

        return cls.datamanager


    @classmethod
    def removeFinishedProcs(cls, vActiveProcs):

        vDelProcs = []

        for i in range(0, len(vActiveProcs)):
            oProc = vActiveProcs[i]

            oProc.join(0.05)

            if not oProc.is_alive():

                if oProc.exitcode != 0:
                    eprint(str(oProc) + " " + str(oProc.exitcode))

                vDelProcs.append(i)

        if len(vDelProcs) > 0:
            vDelProcs = sorted(vDelProcs, reverse=True)

            for x in vDelProcs:
                #eprint("deleting: " + str(x))
                del vActiveProcs[x]

        return vActiveProcs

    @classmethod
    def mapReduce(cls, iThreads, oIterable, oFunc, sEnvironment, chunkSize = 1, pReduceFunc = None):

        """

        Process based parallelism

        :param iThreads:
        :param oIterable:
        :param oFunc:
        :param sEnvironment:
        :param chunkSize:
        :param pReduceFunc:
        :return:
        """

        workQueue = multiprocessing.Queue(len(oIterable))
        resultQueue = multiprocessing.JoinableQueue(len(oIterable))
        oReturnQueueLock = multiprocessing.Lock()
        vReturn = []

        for x in oIterable:
            workQueue.put(x)

        iProcCount = 0

        vReturn = None

        if pReduceFunc is None:
            vReturn = []

        vActiveProcs = []

        bWorkQueueEmpty = workQueue.empty()
        #eprint("WorkQueue size: " + str(workQueue.qsize()) + " and workqueue.empty " + str(bWorkQueueEmpty))

        bWorkQueueEmpty = True
        while bWorkQueueEmpty:

            while len(vActiveProcs) >= iThreads:

                try:

                    while True:

                        iStart = time.time()

                        oElem = resultQueue.get(timeout=0.05)
                        resultQueue.task_done()

                        if not pReduceFunc is None:

                            vReturn = pReduceFunc(vReturn, oElem, sEnvironment)

                            iEnd = time.time()

                            #eprint("Reduced: " + str(iEnd-iStart))
                        else:
                            vReturn.append(oElem)

                except queue.Empty:
                    time.sleep(0.05)

                vActiveProcs = cls.removeFinishedProcs(vActiveProcs)

                continue

            # create chunk
            vChunk = []
            for i in range(0, chunkSize):

                try:
                    oElem = workQueue.get(timeout=0.05)
                    vChunk.append( oElem )
                except queue.Empty:

                    bWorkQueueEmpty = False

                    #eprint("WorkQueue Empty")

                    break


            oProc = Parallel.MyProcess(iProcCount, vChunk, resultQueue, oReturnQueueLock, sEnvironment, oFunc)
            iProcCount += 1

            vActiveProcs.append( oProc )
            oProc.start()
            #eprint("Started " + str(oProc))


        #eprint("Reduce Phase")

        while len(vActiveProcs) > 0:

            try:

                while True:
                    oElem = resultQueue.get(timeout=0.05)
                    resultQueue.task_done()

                    if not pReduceFunc is None:
                        vReturn = pReduceFunc(vReturn, oElem, sEnvironment)
                    else:
                        vReturn.append(oElem)

            except queue.Empty:
                vActiveProcs = cls.removeFinishedProcs(vActiveProcs)

            continue


        return vReturn


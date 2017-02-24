import sys
import os
import glob
import tarfile
import shutil
import h5py
from enum import Enum
from collections import OrderedDict, Counter

import mjoppich.utils.FileUtils as fu


class Fast5TYPE(Enum):
    PRE_BASECALL = 4,
    BASECALL_RNN_1D = 3,
    BASECALL_1D = 2,
    BASECALL_1D_COMPL=1,
    BASECALL_2D=0,

    UNKNOWN=-1

class Fast5PATH(Enum):
    PRE_BASECALL=0,
    TEMPLATE=1,
    COMPLEMENT=2,
    TEMPLATE_2D=3,

class classproperty(object):

    def __init__(self, fget):
        self.fget = fget

    def __get__(self, owner_self, owner_cls):
        return self.fget(owner_cls)

class Fast5File:

    @classproperty
    def paths(cls):
        return OrderedDict([
            (Fast5TYPE.BASECALL_2D,
                {
                    Fast5PATH.TEMPLATE: '/Analyses/Basecall_2D_%03d/BaseCalled_template',
                    Fast5PATH.COMPLEMENT: '/Analyses/Basecall_2D_%03d/BaseCalled_complement',
                    Fast5PATH.TEMPLATE_2D: '/Analyses/Basecall_2D_%03d/BaseCalled_2D',
                    Fast5PATH.PRE_BASECALL: '/Analyses/EventDetection_%03d/Reads/'
                }),
            (Fast5TYPE.BASECALL_1D_COMPL,
             {
                 Fast5PATH.TEMPLATE: '/Analyses/Basecall_1D_%03d/BaseCalled_template',
                 Fast5PATH.COMPLEMENT: '/Analyses/Basecall_1D_%03d/BaseCalled_complement',
                 Fast5PATH.PRE_BASECALL: '/Analyses/EventDetection_%03d/Reads/'
             }),
            (Fast5TYPE.BASECALL_1D,
             {
                 Fast5PATH.TEMPLATE: '/Analyses/Basecall_1D_%03d/BaseCalled_template',
                 Fast5PATH.PRE_BASECALL: '/Analyses/EventDetection_%03d/Reads/'
             }),
            (Fast5TYPE.BASECALL_RNN_1D,
             {
                 Fast5PATH.TEMPLATE: '/Analyses/Basecall_RNN_1D_%03d/BaseCalled_template'
             }),
            (Fast5TYPE.PRE_BASECALL,
             {
                 Fast5PATH.PRE_BASECALL: '/Analyses/EventDetection_%03d/Reads/'
             }),
            (Fast5TYPE.UNKNOWN, {})

        ])

    @classproperty
    def discriminating_paths(cls):
        return OrderedDict([
            (Fast5TYPE.BASECALL_2D,
                {
                    Fast5PATH.TEMPLATE: '/Analyses/Basecall_2D_%03d/BaseCalled_template',
                    Fast5PATH.COMPLEMENT: '/Analyses/Basecall_2D_%03d/BaseCalled_complement',
                    Fast5PATH.TEMPLATE_2D: '/Analyses/Basecall_2D_%03d/BaseCalled_2D',
                }),
            (Fast5TYPE.BASECALL_1D_COMPL,
             {
                 Fast5PATH.COMPLEMENT: '/Analyses/Basecall_1D_%03d/BaseCalled_complement',
             }),
            (Fast5TYPE.BASECALL_1D,
             {
                 Fast5PATH.TEMPLATE: '/Analyses/Basecall_1D_%03d/BaseCalled_template',
             }),
            (Fast5TYPE.BASECALL_RNN_1D,
             {
                 Fast5PATH.TEMPLATE: '/Analyses/Basecall_RNN_1D_%03d/BaseCalled_template'
             }),
            (Fast5TYPE.PRE_BASECALL,
             {
                 Fast5PATH.PRE_BASECALL: '/Analyses/EventDetection_%03d/Reads/'
             }),
            (Fast5TYPE.UNKNOWN, {})

        ])

    def __init__(self, path, group = 0, keep_file_open = False):

        self.group = group
        self.filename = path
        self.type = None

        self.hdf5file = None
        self.is_open = self._open(keep_file_open)

        self.fasta = []
        self.fastq = []

        self.has_fastq = False
        self.has_fasta = False
        self.has_pre_basecalled = False
        self.has_metadata = False

        self.pore = (-1,-1)
        self.timestamp = 0
        self.readnum = 0
        self.read_length = 0

    def _open(self, keep_open):

        self.hdf5file = h5py.File(self.filename, 'r')
        self.type = self._guessType()


    def hdf_error(self, reason):

        msg = """hdf reading failure in file '%s': %s """ % (self.filename, reason)

        sys.exit(msg)

    def _guessType(self):

        for filetype in self.discriminating_paths:

            needed_paths = self.discriminating_paths[filetype]

            vpaths = [needed_paths[x] for x in needed_paths]

            pathsFound = False

            for path in vpaths:

                search_path = path % (self.group)

                if search_path in self.hdf5file:
                    pathsFound = True
                    # finding a single path is sufficient!
                    break

            if pathsFound:
                return filetype

        return Fast5TYPE.UNKNOWN

    def _read_fastq(self, type):

        try:

            avail_paths = self.paths[type]

            avail_fastq = {}

            for path_type in avail_paths:
                path = avail_paths[path_type]

                hdf_path = path % self.group
                hdf_elem = self.hdf5file[hdf_path]

                content = hdf_elem['Fastq'][()]

                avail_fastq[path_type] = content

            return avail_fastq

        except Exception:
            return None


    def _read_attrib(self, attribname):

        try:
            path = "/Raw/Reads/"

            readsGroup = self.hdf5file[path]

            storedReads = readsGroup.keys()

            readNum = -1
            singleRead = True
            if len(storedReads) > 1:
                singleRead = False

                readNum = set()

            for read in storedReads:

                seqRead = self.hdf5file[path + read]

                if singleRead:
                    return seqRead.attrs[ attribname ]
                else:
                    readNum.add(seqRead.attrs[ attribname ])

            return readNum

        except:

            return None

    def readNumber(self):
        """

        :return: read number, or set of read numbers (if multiple reads in file)
        """
        return self._read_attrib('read_number')

    def readID(self):
        """

        :return: read number, or set of read numbers (if multiple reads in file)
        """
        return str(self._read_attrib('read_id'))


    def runID(self):

        try:
            path = "/UniqueGlobalKey/tracking_id/"
            hdf_elem = self.hdf5file[path]

            run_id = hdf_elem.attrs['run_id']

            return str(run_id)

        except:
            return None

        return None

    def channelID(self):

        try:
            path = "/UniqueGlobalKey/channel_id/"
            hdf_elem = self.hdf5file[path]

            channel_number = hdf_elem.attrs['channel_number']

            return str(channel_number)

        except:
            return None

        return None


    def _read_main_fastq(self):

        return self._read_fastq(self.type)

    def analyses(self):

        return ( [str(x) for x in self.hdf5file['Analyses']] )

    def printGroupsAttribs(self):

        def printName(name):



            print(str(name) + " " + str([(k, v) for (k,v) in self.hdf5file[name].attrs.items()]))




        self.hdf5file.visit(printName)



class Fast5Directory:

    def __init__(self, path):

        if not os.path.isdir(path):
            raise ValueError("Given path is not a directory: " + path)

        self.path = fu.makePath(path)
        self.filesIT = glob.iglob(self.path + "*.fast5")

    def collect(self):

        for nextFile in self.filesIT:

            yield Fast5File(nextFile)


if __name__ == "__main__":

    folderpath = "/home/proj/projekte/sequencing/HPYLORI_MINION_ILLUMINA_JIMENEZ/reads/reads_P12_pooled/fail/"

    fileCnt = 0

    for folderpath in ["/home/proj/projekte/sequencing/HPYLORI_MINION_ILLUMINA_JIMENEZ/reads/reads_P12_pooled/fail/",
        "/home/proj/projekte/sequencing/HPYLORI_MINION_ILLUMINA_JIMENEZ/reads/reads_P12_pooled/pass/",
        "/home/proj/projekte/sequencing/HPYLORI_MINION_ILLUMINA_JIMENEZ/reads/reads_P12_pooled/skip/"]:

        fileIt = fu.getFilePatternIterator(folderpath, "*.fast5")

        cntDupReadID = Counter()

        for file in fileIt:
            if "mux_scan" in file:
                continue

            fileCnt += 1

            if fileCnt % 1000 == 0:
                print(fileCnt)

            oFile = Fast5File(file)

            # oFile.printGroupsAttribs()

            cntDupReadID[oFile.readID()] += 1

        for x in cntDupReadID:

            if cntDupReadID[x] > 1:
                print(str(x) + " " + cntDupReadID[x])

    exit(0)

    for x in Fast5File.paths:
        print(x)
        print(Fast5File.paths[x])

    folderpath = "/home/proj/projekte/sequencing/HPYLORI_MINION_ILLUMINA_JIMENEZ/reads/fail/uploaded/"

    fileIt = fu.getFilePatternIterator(folderpath, "MVPI22502_20161222_FNFAB31358_MN19136_sequencing_run_2d*.fast5")

    countSuccess = Counter()
    countAnalyses = Counter()
    countBases2D = Counter()

    totalFiles = 0
    fastqCounts = 0

    for file in fileIt:

        totalFiles += 1

        try:
            oFile = Fast5File(file)

            countSuccess[ oFile.type ] += 1
        except:

            countSuccess['exception'] += 1
            continue

        if oFile.type == Fast5TYPE.UNKNOWN:
            print(file)

        for analysis in oFile.analyses():
            countAnalyses [analysis] += 1

        if 'Basecall_2D_000' in oFile.analyses():

            for x in oFile.hdf5file['/Analyses/Basecall_2D_000/']:
                countBases2D[x]+=1

            if "/Analyses/Basecall_2D_000/BaseCalled_2D" in oFile.hdf5file:

                fastq = oFile.hdf5file['/Analyses/Basecall_2D_000/BaseCalled_2D/']['Fastq'][()]

                if fastq != None:
                    #print( fastq.decode('UTF-8') )
                    fastqCounts += 1

            #print(oFile._read_fastq(Fast5TYPE.BASECALL_2D))


    print(countSuccess)
    print(countAnalyses)
    print(countBases2D)
    print(totalFiles)
    print(fastqCounts)


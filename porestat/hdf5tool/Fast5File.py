import sys
import os
import glob
import tarfile
import shutil
import h5py
from enum import Enum
from collections import OrderedDict, Counter

import mjoppich.utils.FileUtils as fu

from porestat.hdf5tool.SequenceFormats import FASTQ

class classproperty(object):

    def __init__(self, fget):
        self.fget = fget

    def __get__(self, owner_self, owner_cls):
        return self.fget(owner_cls)

class Fast5TYPE(Enum):
    PRE_BASECALL = 4,
    BASECALL_RNN_1D = 3,
    BASECALL_1D = 2,
    BASECALL_1D_COMPL=1,
    BASECALL_2D=0,

    UNKNOWN=-1

    @classproperty
    def type2str(cls):
        return {

            Fast5TYPE.BASECALL_2D: 'BASECALL_2D',
            Fast5TYPE.BASECALL_1D_COMPL: 'BASECALL_1D_COMPL',
            Fast5TYPE.BASECALL_1D: 'BASECALL_1D',
            Fast5TYPE.BASECALL_RNN_1D: 'BASECALL_RNN_1D',
            Fast5TYPE.PRE_BASECALL: 'PRE_BASECALL',
            Fast5TYPE.UNKNOWN: 'UNKNOWN'
        }

    @classproperty
    def str2type(cls):
        return {cls.type2str[x]: x for x in cls.type2str}

class Fast5PATH(Enum):
    PRE_BASECALL=0,
    TEMPLATE=1,
    COMPLEMENT=2,
    TEMPLATE_2D=3,



class Fast5File:

    @classmethod
    def join_paths(cls, prefix, suffix):

        if prefix[len(prefix)-1] == '/':
            prefix = prefix[0:len(prefix)-1]

        if len(suffix) > 0 and suffix[0] == '/':
            suffix = suffix[1:]

        return prefix + "/" + suffix


    @classproperty
    def analyses_paths(cls):
        return OrderedDict([
            (Fast5TYPE.BASECALL_2D, '/Analyses/Basecall_2D_%03d/'),
            (Fast5TYPE.BASECALL_1D_COMPL, '/Analyses/Basecall_1D_%03d/'),
            (Fast5TYPE.BASECALL_1D, '/Analyses/Basecall_1D_%03d/'),
            (Fast5TYPE.BASECALL_RNN_1D, '/Analyses/Basecall_RNN_1D_%03d/'),
            (Fast5TYPE.PRE_BASECALL, '/Analyses/EventDetection_%03d/')
            #(Fast5TYPE.UNKNOWN, "")
        ])

    @classproperty
    def sequence_paths(cls):
        return OrderedDict([

            (Fast5TYPE.BASECALL_2D, 'BaseCalled_2D'),
            (Fast5TYPE.BASECALL_1D_COMPL, 'BaseCalled_complement'),
            (Fast5TYPE.BASECALL_1D, 'BaseCalled_template'),
            (Fast5TYPE.BASECALL_RNN_1D, 'BaseCalled_template'),
            (Fast5TYPE.PRE_BASECALL, ''),
            (Fast5TYPE.UNKNOWN, "")

        ])

    def __init__(self, path, group = 0, keep_file_open = False):

        self.filename = path
        self.type = None

        self.sequence_paths = {}
        self.winner = Fast5PATH.TEMPLATE

        self.hdf5file = None
        self.is_open = self._open(keep_file_open)

        self.pore = (-1,-1)
        self.timestamp = 0
        self.readnum = 0
        self.read_length = 0

    def _open(self, keep_open):

        self.hdf5file = h5py.File(self.filename, 'r')
        self.type = self._guessType()

        #self.printGroupsAttribs()



    def hdf_error(self, reason):

        msg = """hdf reading failure in file '%s': %s """ % (self.filename, reason)

        sys.exit(msg)

    def _guessType(self):

        for filetype in Fast5File.analyses_paths:

            analyses_path = Fast5File.analyses_paths[filetype]

            group = -1


            # findest highest group for filetype
            while True:

                analyses_group = analyses_path % (group+1)

                if not analyses_group in self.hdf5file:
                    break

                group += 1

            if group == -1:
                continue

            analyses_group = analyses_path % group
            seqpath = self.join_paths(analyses_group, Fast5File.sequence_paths[filetype])

            if not seqpath in self.hdf5file:
                continue

            self.sequence_paths = {}

            # get analyses_group . 'basecall_1d'
            if filetype == Fast5TYPE.BASECALL_2D:

                basecall1dpath = '/' + self._get_attribute(analyses_group, 'basecall_1d')

                self.sequence_paths[Fast5PATH.TEMPLATE] = self.join_paths( basecall1dpath , "BaseCalled_template" )
                self.sequence_paths[Fast5PATH.COMPLEMENT] = self.join_paths( basecall1dpath , "BaseCalled_complement" )
                self.sequence_paths[Fast5PATH.TEMPLATE_2D] = seqpath

                eventpath = '/' + self._get_attribute(basecall1dpath, 'event_detection', "Analyses/EventDetection_000")
                self.sequence_paths[Fast5PATH.PRE_BASECALL] = eventpath

                self.winner = Fast5PATH.TEMPLATE_2D

            elif filetype == Fast5TYPE.BASECALL_1D_COMPL:

                self.sequence_paths[Fast5PATH.TEMPLATE] = self.join_paths(analyses_group, "BaseCalled_template")
                self.sequence_paths[Fast5PATH.COMPLEMENT] = self.join_paths(analyses_group, "BaseCalled_complement")

                eventpath = '/' + self._get_attribute(analyses_group, 'event_detection',"Analyses/EventDetection_000")
                self.sequence_paths[Fast5PATH.PRE_BASECALL] = eventpath

                self.winner = Fast5PATH.COMPLEMENT

            elif filetype == Fast5TYPE.BASECALL_1D:

                self.sequence_paths[Fast5PATH.TEMPLATE] = self.join_paths(analyses_group, "BaseCalled_template")

                eventpath = '/' + self._get_attribute(analyses_group, 'event_detection', "Analyses/EventDetection_000")
                self.sequence_paths[Fast5PATH.PRE_BASECALL] = eventpath

                self.winner = Fast5PATH.TEMPLATE


            elif filetype == Fast5TYPE.BASECALL_RNN_1D:

                self.sequence_paths[Fast5PATH.TEMPLATE] = self.join_paths(analyses_group, "BaseCalled_template")

                eventpath = '/' + self._get_attribute(analyses_group, 'event_detection', "Analyses/EventDetection_000")
                self.sequence_paths[Fast5PATH.PRE_BASECALL] = eventpath

                self.winner = Fast5PATH.TEMPLATE


            elif filetype == Fast5TYPE.PRE_BASECALL:

                self.sequence_paths[Fast5PATH.PRE_BASECALL] = seqpath

                self.winner = None

            # set up available sequences

            # sanity check
            for path in self.sequence_paths:

                temp_path = self.sequence_paths[path]
                exists = temp_path in self.hdf5file

                if not exists:
                    print( temp_path + " " + str(temp_path in self.hdf5file))

            return filetype

        return Fast5TYPE.UNKNOWN

    def getFastQ(self, type = None):

        if type == None:
            type = self.winner

        return self._read_fastq(type)

    def getFastA(self, type = None):

        if type == None:
            type = self.winner

        fastq = self._read_fastq(type)

        if fastq != None:
            return fastq.to_fasta()

        return None

    def _read_fastq(self, type):

        try:

            if type == None:
                return None

            pathToFQ = self.sequence_paths[type]

            hdf_elem = self.hdf5file[pathToFQ]
            fastq_seq = hdf_elem['Fastq'][()].decode("utf-8")

            return FASTQ.parseFromStr(fastq_seq)

        except Exception as e:

            if self.runID() != "84f1d5bc960cc73afb857b227c377da6f10ce650":
                return None

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

    def _get_attribute(self, path, attrib, default = None):
        try:
            hdf_elem = self.hdf5file[path]

            run_id = hdf_elem.attrs[attrib]

            return run_id.decode("utf-8")

        except:
            return default

    def user_filename_input(self):
        """

        :return: user_filename_input from UniqueGlobalKey/context_tags
        """
        return self._get_attribute("/UniqueGlobalKey/context_tags/", 'user_filename_input')

    def readNumber(self):
        """

        :return: read number, or set of read numbers (if multiple reads in file)
        """
        return int(self._read_attrib('read_number'))

    def readID(self):
        """

        :return: read number, or set of read numbers (if multiple reads in file)
        """
        return self._read_attrib('read_id').decode("utf-8")


    def runID(self):

        return self._get_attribute("/UniqueGlobalKey/tracking_id/", 'run_id')

    def channelID(self):

        return int(self._get_attribute("/UniqueGlobalKey/channel_id/", 'channel_number'))


    def _read_main_fastq(self):

        return self._read_fastq(self.type)

    def analyses(self):

        return ( [str(x) for x in self.hdf5file['Analyses']] )

    def printGroupsAttribs(self):

        def printName(name):
            print(str(name) + " " + str([(k, v) for (k,v) in self.hdf5file[name].attrs.items()]))



        print(self.filename)
        print(self.type)
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


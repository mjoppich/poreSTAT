import argparse
import sys
import os
import glob

import h5py
from enum import Enum
from collections import OrderedDict, Counter


from .SequenceFormats import FASTQ
import dateutil.parser

class classproperty(object):
    """
    class to allow the usage of class properties (extension of class functions)
    """

    def __init__(self, fget):
        self.fget = fget

    def __get__(self, owner_self, owner_cls):
        return self.fget(owner_cls)


class Fast5TYPEAction(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):

        try:

            if len(values) == 0:
                style = None

            else:

                style = []
                if type(values) == list:

                    for elem in values:
                        style.append(Fast5TYPE[elem])
                else:
                    style = [Fast5TYPE[values]]

            args.__dict__[self.dest] = style

        except:

            raise argparse.ArgumentError(None, 'Fast5TYPE can not be {n}, it must be one of {m}'.format(n=values, m=', '.join([str(x.value) for x in Fast5TYPE])))

class Fast5TYPE(Enum):
    PRE_BASECALL = 'PRE_BASECALL'
    BASECALL_RNN_1D = 'BASECALL_RNN_1D'
    BASECALL_1D = 'BASECALL_1D'
    BASECALL_1D_COMPL = 'BASECALL_1D_COMPL'
    BASECALL_2D = 'BASECALL_2D'
    BARCODING = 'BARCODING'
    UNKNOWN = 'UNKNOWN'

    def __str__(self):
        return self.name

    # @classproperty
    # def type2str(cls):
    #     """
    #     :return :dictionary of type to str for fast5 read types
    #     """
    #     return {x: x.name for x in list(cls)}
    #
    # @classproperty
    # def str2type(cls):
    #     """
    #     :return : dictionary of str to type for fast5 read types
    #     """
    #     return {x.name: x for x in list(cls)}

class Fast5FileException(Exception):
    pass
class Fast5FileTYPE(Enum):
    SINGLE = 'SINGLE_FASTA'
    MULTIFASTA = 'MULTI_FASTA'

    def __str__(self):
        return self.name

class MFast5File:

    def _open(self, keep_open):

        self.hdf5file = h5py.File(self.filename, 'r')
        self.type = self._guessType()

    def __init__(self, path, keep_file_open = False):

        self.hdf5file = None
        self.filename = path
        self.type = None
        self.iterPaths = []
        self.is_open = self._open(keep_file_open)

        #print(self.type)

    def _guessType(self):

        cntGrpsMulti = 0
        cntGrpsSingle = 0

        allIterPaths = list()
        for grp in self.hdf5file:

            if grp.startswith("read_"):
                cntGrpsMulti += 1
                if not grp in allIterPaths:
                    allIterPaths.append(grp)
            else:
                cntGrpsSingle += 1

        if cntGrpsSingle == 0:
            self.iterPaths = list(allIterPaths)
            return Fast5FileTYPE.MULTIFASTA

        self.iterPaths = ["/"]
        return Fast5FileTYPE.SINGLE

    def __len__(self):
        return len(self.iterPaths)
        
    def __iter__(self):

        self.iterIdx = 0
        return self

    def __next__(self):

        if self.iterIdx >= len(self.iterPaths):
            raise StopIteration

        iPath = self.iterPaths[self.iterIdx]
        self.iterIdx += 1
        return Fast5File(self.hdf5file[iPath], path=iPath)


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
            (Fast5TYPE.BASECALL_2D, './Analyses/Basecall_2D_%03d/'),
            (Fast5TYPE.BASECALL_1D_COMPL, './Analyses/Basecall_1D_%03d/'),
            (Fast5TYPE.BASECALL_1D, './Analyses/Basecall_1D_%03d/'),
            (Fast5TYPE.BASECALL_RNN_1D, './Analyses/Basecall_RNN_1D_%03d/'),
            (Fast5TYPE.BARCODING, './Analyses/Barcoding_%03d/'),
            (Fast5TYPE.PRE_BASECALL, './Analyses/EventDetection_%03d/')
            #(Fast5TYPE.UNKNOWN, "")
        ])

    @classproperty
    def sequence_paths(cls):
        return OrderedDict([

            (Fast5TYPE.BASECALL_2D, 'BaseCalled_2D'),
            (Fast5TYPE.BASECALL_1D_COMPL, 'BaseCalled_complement'),
            (Fast5TYPE.BASECALL_1D, 'BaseCalled_template'),
            (Fast5TYPE.BASECALL_RNN_1D, 'BaseCalled_template'),
            (Fast5TYPE.BARCODING, 'Barcoding'),
            (Fast5TYPE.PRE_BASECALL, ''),
            (Fast5TYPE.UNKNOWN, "")

        ])

    def __init__(self, h5grp, path=None, group=0, keep_file_open = False):

        self.type = None
        self.filename = path

        self.sequence_paths = {}
        self.winner = Fast5TYPE.BASECALL_1D

        self.hdf5file = h5grp
        self.type = self._guessType()

        self.pore = (-1, -1)
        self.timestamp = 0
        self.readnum = 0
        self.read_length = 0


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
                #try alternative version
                seqpath = self.join_paths(analyses_group, "BaseCalled_template")

                if not seqpath in self.hdf5file:
                    continue

            self.sequence_paths = {}

            # get analyses_group . 'basecall_1d'
            if filetype == Fast5TYPE.BASECALL_2D:

                
                basecallAttrib = self._get_attribute(analyses_group, 'basecall_1d')                
                # if that returns None, attrib not existant => no 2D basecall
                if basecallAttrib != None:
                    
                    basecall1dpath = basecallAttrib #basecall1dpath = '/Analyses/' + basecallAttrib
                    self.sequence_paths[Fast5TYPE.BASECALL_1D] = self.join_paths( basecall1dpath , "BaseCalled_template" )
                    self.sequence_paths[Fast5TYPE.BASECALL_1D_COMPL] = self.join_paths( basecall1dpath , "BaseCalled_complement" )

                    eventpath = '/' + self._get_attribute(basecall1dpath, 'event_detection', "Analyses/EventDetection_000")
                    self.sequence_paths[Fast5TYPE.PRE_BASECALL] = eventpath


                self.sequence_paths[Fast5TYPE.BASECALL_2D] = seqpath
                self.winner = Fast5TYPE.BASECALL_2D

            elif filetype == Fast5TYPE.BASECALL_1D_COMPL:

                self.sequence_paths[Fast5TYPE.BASECALL_1D] = self.join_paths(analyses_group, "BaseCalled_template")
                self.sequence_paths[Fast5TYPE.BASECALL_1D_COMPL] = self.join_paths(analyses_group, "BaseCalled_complement")

                eventpath = '/' + self._get_attribute(analyses_group, 'event_detection',"Analyses/EventDetection_000")
                self.sequence_paths[Fast5TYPE.PRE_BASECALL] = eventpath

                self.winner = Fast5TYPE.BASECALL_1D_COMPL

            elif filetype == Fast5TYPE.BASECALL_1D:

                self.sequence_paths[Fast5TYPE.BASECALL_1D] = self.join_paths(analyses_group, "BaseCalled_template")

                eventpath = '/' + self._get_attribute(analyses_group, 'event_detection', "Analyses/EventDetection_000")
                self.sequence_paths[Fast5TYPE.PRE_BASECALL] = eventpath

                self.winner = Fast5TYPE.BASECALL_1D


            elif filetype == Fast5TYPE.BASECALL_RNN_1D:

                self.sequence_paths[Fast5TYPE.BASECALL_RNN_1D] = self.join_paths(analyses_group, "BaseCalled_template")

                eventpath = '/' + self._get_attribute(analyses_group, 'event_detection', "Analyses/EventDetection_000")
                self.sequence_paths[Fast5TYPE.PRE_BASECALL] = eventpath
                self.winner = Fast5TYPE.BASECALL_RNN_1D


            elif filetype == Fast5TYPE.PRE_BASECALL:

                self.sequence_paths[Fast5TYPE.PRE_BASECALL] = seqpath
                self.winner = None

            # set up available sequences

            # sanity check
            ftypeEjected=False
            for path in self.sequence_paths:

                temp_path = self.sequence_paths[path]
                exists = temp_path in self.hdf5file

                if not exists and not path == Fast5TYPE.PRE_BASECALL:
                    #print( temp_path + " " + str(temp_path in self.hdf5file))
                    ftypeEjected=True
                    
            if ftypeEjected:
                #print("Eject ftype", filetype)
                continue

            return filetype

        return Fast5TYPE.UNKNOWN

    def getFastQ(self, readType=None):

        if readType == None:
            readType = self.winner
       
        return self._read_fastq(readType)


    def _get_signal(self):

        iReadNum = self.readNumber()

        sSignalPath = 'Raw/Reads/Read_' + str(iReadNum)

        hdf_elem = self.hdf5file[sSignalPath]

        raw_signal = hdf_elem['Signal'][()]

        iSignalLength = len(raw_signal)

        return iSignalLength

    def _read_fastq(self, readType):

        try:

            if readType == None:
                return None

            pathToFQ = self.sequence_paths[readType]

            hdf_elem = self.hdf5file[pathToFQ]
            fastq_seq = hdf_elem['Fastq'][()]

            if not type(fastq_seq) == str:
                fastq_seq = fastq_seq.decode("utf-8")

            return FASTQ.parseFromStr(fastq_seq)

        except Exception as e:

            return None

    def _read_raw_signal(self):

        try:
            allevs=[]
            readNames = [x for x in self.hdf5file['Raw/Reads']]

            for readName in readNames:
                hdf_elem = self.hdf5file['Raw/Reads/'+readName]
                evs = hdf_elem['Signal'][()]

                allevs.append(evs)

            return allevs



        except Exception as e:

            return []

    def _read_events(self, readType):

        try:

            if readType == None:
                return None

            pathToFQ = self.sequence_paths[readType]

            hdf_elem = self.hdf5file[pathToFQ]
            fastq_seq = hdf_elem['Events'][()]

            return fastq_seq

        except Exception as e:

            return None

    def getSampleFrequency(self):
        
        if "context_tags" in self.hdf5file:
            return int(self._get_attribute('context_tags', 'sample_frequency', 4000)) 
        else:
            return int(self._get_attribute('/UniqueGlobalKey/context_tags', 'sample_frequency', 4000))

    def getExperimentStartTime(self):


        if "tracking_id" in self.hdf5file:
            timeAttribData = self._get_attribute('tracking_id', 'exp_start_time', None)
        else:
            timeAttribData = self._get_attribute('/UniqueGlobalKey/tracking_id', 'exp_start_time', None)

        try:
            timestamp = int(timeAttribData)
            return timestamp
        except:
            try:
                
                parsedTime = dateutil.parser.parse(str(timeAttribData))
                timestamp = parsedTime.timestamp()
                return int(timestamp)

            except:
                raise Fast5FileException("Invalid exp start time format: " + str(timeAttribData))



    def readCreateTime(self):


        expStartTime = self.getExperimentStartTime()
        sampleFrequency = self.getSampleFrequency()

        #print(expStartTime)
        #print(sampleFrequency)

        if expStartTime == None:
            return None

        #self.printGroupsAttribs()
        startTimes = self._read_attrib('start_time')

        #print(startTimes)

        if startTimes == None:
            path = "/Analyses/EventDetection_000/Reads/"+self.readID()
            startTimes = self._get_attribute(path, "start_time", None)
            

        if startTimes == None:
            return None

        if type(startTimes) == list:
            vReturn = [x + expStartTime for x in startTimes]
            return vReturn

        else:
            iReadStart = int(expStartTime + (startTimes / sampleFrequency))
            return iReadStart


    def _read_attrib(self, attribname):

        try:
            path = "/Raw/Reads/"

            if not path in self.hdf5file:
                path = "/Analyses/EventDetection_000/Reads/"    

            if not path in self.hdf5file:
                path = "Raw"    


            if path != "Raw":

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
                        return seqRead.attrs[attribname]
                    else:
                        readNum.add(seqRead.attrs[attribname])

                return readNum
            
            else:
                
                readsGroup = self.hdf5file[path]
                return readsGroup.attrs[attribname]


        except:

            return None

    def _get_attribute(self, path, attrib, default=None):
        try:

            if path in self.hdf5file:
                hdf_elem = self.hdf5file[path]

                run_id = hdf_elem.attrs[attrib]

                if not type(run_id) == str:
                    return run_id.decode("utf-8")
                else:
                    return run_id

            raise Fast5FileException("Path does not exist: " + path)

        except:
            return default

    def user_filename_input(self):
        """
        :return: user_filename_input from UniqueGlobalKey/context_tags
        """
        if "context_tags" in self.hdf5file:
            if "user_filename_input" in self.hdf5file["context_tags"].attrs:
                return self._get_attribute('context_tags', 'user_filename_input', None)
            else:
                return self._get_attribute('tracking_id', 'sample_id', None)
        else:

            if "user_filename_input" in self.hdf5file["/UniqueGlobalKey/context_tags"].attrs:
                return self._get_attribute('/UniqueGlobalKey/context_tags', 'user_filename_input', None)
            else:
                return self._get_attribute('/UniqueGlobalKey/tracking_id', 'sample_id', None)
            

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
        if "tracking_id" in self.hdf5file:
            return self._get_attribute('tracking_id', 'run_id')
        else:
            return self._get_attribute('/UniqueGlobalKey/tracking_id', 'run_id')


    def channelID(self):

        if "channel_id" in self.hdf5file:
            return int(self._get_attribute('channel_id', 'channel_number', 4000)) 
        else:
            return int(self._get_attribute('/UniqueGlobalKey/channel_id', 'channel_number', 4000))


    def _read_main_fastq(self):

        return self._read_fastq(self.type)

    def analyses(self):

        return [str(x) for x in self.hdf5file['Analyses']]

    def printGroupsAttribs(self):

        def printName(name):
            print(str(name) + " " + str([(k, v) for (k, v) in self.hdf5file[name].attrs.items()]))

        print(self.filename)
        print(self.type)

        self.hdf5file.visit(printName)

#>tig00000001 len=4527241 reads=33593 covStat=14412.82 gappedBases=no class=contig suggestRepeat=no suggestCircular=no

    def sequenceLength(self):

        seq = self.getFastQ()

        if seq == None:
            return -1

        return len(seq)

    def sequenceName(self):

        seq = self.getFastQ()

        if seq == None:
            return ""

        return seq.id

class Fast5Directory:

    def __init__(self, path):

        if not os.path.isdir(path) and not os.path.exists(path):
            raise ValueError("Given path is not a directory or does not exist: " + path)

        self.path = os.path.abspath(path)#Path(path)
        
        self.filesIT = None
        if os.path.isdir(path):
            self.filesIT = glob.glob(self.path +"/" + "**/*.fast5", recursive=True)
            self.filesIT.sort(key=os.path.getctime)
        else:
            self.filesIT = [self.path]

    def collect(self):

        for nextFile in self.filesIT:

            f5file = MFast5File(nextFile)

            if f5file.type == Fast5FileTYPE.MULTIFASTA:
                print("Loaded", nextFile)

            for read in f5file:
                yield read

if __name__ == "__main__":

    testFile = "/home/proj/projekte/sequencing/HPYLORI_MINION_ILLUMINA_JIMENEZ/reads/minion_p12_030417/chip2/reads_170425_chip2_run1_Tx30a_bc/workspace/0/"
    testFile += "MVPI22502_20170425_FNFAF18023_MN19136_sequencing_run_20170425_Sequencing_run_2_1_Tx30a_pooled_57075_ch24_read133_strand.fast5"
    f5File = Fast5File(testFile)



    exit(0)

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


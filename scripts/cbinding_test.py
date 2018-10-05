import ctypes
from collections import Counter

lib = ctypes.cdll.LoadLibrary('../clib/lib/libfoo.so')

class COUNTER_PAIR(ctypes.Structure):
    _fields_ = [
        ("key", ctypes.c_char_p),
        ("value", ctypes.c_uint32)
    ]

class SEQ_COVERAGE(ctypes.Structure):
    _fields_ = [
        ("seq", ctypes.c_char_p),
        ("start", ctypes.c_uint32),
        ("end", ctypes.c_uint32)
    ]

class NT_SUBSTITUTION(ctypes.Structure):
    _fields_ = [
        ("src", ctypes.c_char),
        ("tgt", ctypes.c_char),
        ("count", ctypes.c_uint32)
    ]

class READ_ALIGNMENT(ctypes.Structure):
    _fields_ = [
        ("pReadID", ctypes.c_char_p),
        ("iReadLength", ctypes.c_uint32),
        ("iRefLength", ctypes.c_uint32),
        ("iAlignQual", ctypes.c_uint32),
        ("fSeqQual_min", ctypes.c_float),
        ("fSeqQual_median", ctypes.c_float),
        ("fSeqQual_max", ctypes.c_float),
        ("iLongestMatched", ctypes.c_uint32),
        ("fReadGCContent", ctypes.c_float),
        ("fRefGCContent", ctypes.c_float),
        ("fReadIdentity", ctypes.c_float),
        ("fRefIdentity", ctypes.c_float),

        ("PERF_KMER_COUNT", ctypes.c_uint32),
        ("PERF_KMERS", ctypes.c_void_p),

        ("CIGAR_COUNT", ctypes.c_uint32),
        ("CIGARS", ctypes.c_char_p),
        ("CIGAR_VEC_LENGTHS", ctypes.POINTER(ctypes.c_uint32)),
        ("CIGAR_LENGTHS", ctypes.POINTER(ctypes.POINTER(ctypes.c_uint32))),

        ("MM2KMER_VEC_LENGTHS", ctypes.POINTER(ctypes.c_uint32)),
        ("MM2KMER_LENGTHS", ctypes.POINTER(ctypes.POINTER(ctypes.c_char_p))),

        ("COV_COUNT", ctypes.c_uint32),
        ("COVERAGES", ctypes.c_void_p),

        ("NT_SUBST_COUNT", ctypes.c_uint32),
        ("NT_SUBST", ctypes.c_void_p),

    ]


class AlignmentStatistics(object):
    def __init__(self):
        lib.AlignmentStatistics_new.argtypes = []
        lib.AlignmentStatistics_new.restype = ctypes.c_void_p

        lib.AlignmentStatistics_process.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
        lib.AlignmentStatistics_process.restype = ctypes.c_void_p

        lib.AlignmentStatistics_load_fasta.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        lib.AlignmentStatistics_load_fasta.restype = ctypes.c_void_p

        lib.AlignmentStatistics_ReadStats.argtypes = [ctypes.c_void_p]
        lib.AlignmentStatistics_ReadStats.restype = ctypes.c_void_p

        lib.AlignmentStatistics_ReadStats_size.argtypes = [ctypes.c_void_p]
        lib.AlignmentStatistics_ReadStats_size.restype = ctypes.c_uint32

        self.obj = lib.AlignmentStatistics_new()

    def bar(self):
        lib.Foo_bar(self.obj)
        #
    def loadFasta(self, fastaPath):

        fp = fastaPath.encode('utf-8')

        charStr = ctypes.cast(fp, ctypes.c_char_p)

        lib.AlignmentStatistics_load_fasta(self.obj, charStr)

    def processFiles(self, vals):


        vals = [x.encode('utf-8') for x in vals]
        c_array = (ctypes.c_char_p * len(vals))(*vals)

        charPArray = ctypes.cast(c_array, ctypes.POINTER(ctypes.c_char_p))

        print(type(charPArray))

        retSize = lib.AlignmentStatistics_process(self.obj, ctypes.c_int(len(vals)), charPArray)
        print(retSize)

        retVecSize = lib.AlignmentStatistics_ReadStats_size(self.obj)
        print("Vector Size")
        print(retVecSize)


        retVec = lib.AlignmentStatistics_ReadStats(self.obj)

        s = ctypes.cast(retVec, ctypes.POINTER(READ_ALIGNMENT*retVecSize))

        allElems = [x for x in s.contents]

        print("Printing vector elems")
        for i in range(0, retVecSize):

            elem = allElems[i]

            """
            
            CIGAR2LEN
            """
            cigars = [x for x in elem.CIGARS.decode()[0:elem.CIGAR_COUNT]]
            vecLengths = [x for x in ctypes.cast(elem.CIGAR_VEC_LENGTHS, ctypes.POINTER(ctypes.c_uint32*elem.CIGAR_COUNT)).contents]
            cigarVecs = [x for x in ctypes.cast(elem.CIGAR_LENGTHS, ctypes.POINTER(ctypes.c_void_p*elem.CIGAR_COUNT)).contents]


            cigar2len = {}
            for c in range(0, len(cigars)):

                cigarLen = [x for x in ctypes.cast(cigarVecs[c], ctypes.POINTER(ctypes.c_uint32*vecLengths[c])).contents]

                cigar2len[cigars[c]] = cigarLen

            print("CIGAR2LEN")
            for x in cigar2len:
                print(x, cigar2len[x])


            """
            
            PERF KMERS
            """
            perfKmers = [x for x in
                         ctypes.cast(elem.PERF_KMERS, ctypes.POINTER(COUNTER_PAIR * elem.PERF_KMER_COUNT)).contents]

            perfKmerCounter = Counter()
            for x in perfKmers:
                perfKmerCounter[x.key.decode()] = x.value

            print("PERF KMER")
            for x in perfKmerCounter:
                print(x, perfKmerCounter[x])

            """

            CIGAR2LEN
            """
            mmLengths = [x for x in ctypes.cast(elem.MM2KMER_VEC_LENGTHS,
                                                 ctypes.POINTER(ctypes.c_uint32 * elem.CIGAR_COUNT)).contents]
            mm2kmers = [x for x in
                         ctypes.cast(elem.MM2KMER_LENGTHS, ctypes.POINTER(ctypes.c_void_p * elem.CIGAR_COUNT)).contents]

            mm2kmer = {}
            for c in range(0, len(cigars)):
                mmKmers = [x.decode() for x in
                            ctypes.cast(mm2kmers[c], ctypes.POINTER(ctypes.c_char_p * mmLengths[c])).contents]

                mm2kmer[cigars[c]] = mmKmers

            print("MM2KMER")
            for x in mm2kmer:
                print(x, mm2kmer[x])

            """

            COVERAGES
            """
            covIntervals = [x for x in
                         ctypes.cast(elem.COVERAGES, ctypes.POINTER(SEQ_COVERAGE * elem.COV_COUNT)).contents]

            # TODO add this to global array

            """

            NT SUBSTITUTION
            """
            ntSubst = [x for x in
                         ctypes.cast(elem.NT_SUBST, ctypes.POINTER(NT_SUBSTITUTION * elem.NT_SUBST_COUNT)).contents]

            substCounter = Counter()
            for x in ntSubst:
                substCounter[(x.src.decode(), x.tgt.decode())] = x.count


            print("NT SUBST")
            for x in substCounter:
                print(x, substCounter[x])

            print(allElems[i].pReadID)

        return None




if __name__ == '__main__':
    f = AlignmentStatistics()

    print("Loading FASTA")
    f.loadFasta("/mnt/c/igem/NES/NES.phage.fasta")

    print("Processing SAM")
    print(f.processFiles(['/mnt/c/igem/NES/NES.phage.small.sam']))
import ctypes

lib = ctypes.cdll.LoadLibrary('../clib/build/libfoo.so')


class AlignmentStatistics(object):
    def __init__(self):
        lib.AlignmentStatistics_new.argtypes = []
        lib.AlignmentStatistics_new.restype = ctypes.c_void_p

        lib.AlignmentStatistics_process.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
        lib.AlignmentStatistics_process.restype = ctypes.c_void_p

        self.obj = lib.AlignmentStatistics_new()

    def bar(self):
        lib.Foo_bar(self.obj)

    def processFiles(self, vals):


        vals = [x.encode('utf-8') for x in vals]
        c_array = (ctypes.c_char_p * len(vals))(*vals)

        charPArray = ctypes.cast(c_array, ctypes.POINTER(ctypes.c_char_p))

        print(type(charPArray))

        return lib.AlignmentStatistics_process(self.obj, ctypes.c_int(len(vals)), charPArray)




if __name__ == '__main__':
    f = AlignmentStatistics()

    print(f.processFiles(['/mnt/c/igem/20180810_0938_201808010_NES_500ng/alignment.sam']))
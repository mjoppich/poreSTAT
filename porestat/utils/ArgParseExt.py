import argparse
import os
import tempfile
import shutil
import atexit

from argparse import ArgumentTypeError as err
import os

class FileStubType(object):
    def __init__(self, mode):

        assert mode in ('r','w','t')# or hasattr(type,'__call__')

        self._mode = mode

    def __call__(self, values):

        prospective_dir=values

        prospective_path = os.path.dirname(prospective_dir)

        if 'r' in self._mode:

            if not os.path.isdir(prospective_path):
                raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
            if os.access(prospective_path, os.R_OK):
                return prospective_dir
            else:
                raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

        elif 'w' in self._mode:

            if os.path.isdir(prospective_path):
                if os.access(prospective_path, os.W_OK):
                    return prospective_dir
                else:
                    raise argparse.ArgumentTypeError("writable_dir:{0} is not a writable dir".format(prospective_dir))



class FolderType(object):
    def __init__(self, mode):

        assert mode in ('r','w','t')# or hasattr(type,'__call__')

        self._mode = mode

    def __call__(self, values):

        prospective_dir=values

        if 'r' in self._mode:

            if not os.path.isdir(prospective_dir):
                raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
            if os.access(prospective_dir, os.R_OK):
                return prospective_dir
            else:
                raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

        elif 'w' in self._mode:

            if os.path.isdir(prospective_dir):
                if os.access(prospective_dir, os.W_OK):
                    return prospective_dir
                else:
                    raise argparse.ArgumentTypeError("readable_dir:{0} is not a writable dir".format(prospective_dir))

            else:
                os.makedirs(prospective_dir, exist_ok=True)

                return prospective_dir

        raise argparse.ArgumentTypeError("Directory Input:{0} is not a valid dir".format(prospective_dir))



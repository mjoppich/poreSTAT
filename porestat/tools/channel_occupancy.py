from .PTToolInterface import PTToolInterface
from ..hdf5tool.Fast5File import Fast5File

class Channel_occupancy(PTToolInterface):

    def __init__(self, parser, subparsers):

        super(Channel_occupancy, self).__init__(parser, self.__addParser(subparsers))


    def __addParser(self, subparsers):

        parser_chocc = subparsers.add_parser('occ', help='occ help')
        parser_chocc.add_argument('folders', nargs='+', type=str, help='folders to scan')
        parser_chocc.add_argument('-r', '--reads', type=str, help='minion read folder', required=False)
        parser_chocc.set_defaults(func=self.exec)


        return parser_chocc


    def exec(self, args):

        print("executed CO")